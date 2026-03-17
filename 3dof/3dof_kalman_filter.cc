#include "3dof_kalman_filter.h"
#include "3dof_kalman_filter_mat.h"
#include <cmath>
#include <cstdlib>
#include <cstring>

namespace mobili::ins {

static filter* filter_;
static FilterParams g_params;

// Gravity [m/s^2]
static const float kG[3] = {0.f, 0.f, 9.81f};

void ComputeObservation(const float R_wc[9], const float R_cm[9],
                        float dv_dt, float velocity_cur,
                        const float gyr_imu[3], float h[3]) {
  float* RwcT = filter_->tmp_3x3_3;
  float* RcmT = filter_->tmp_3x3_4;
  float* RwcTg = filter_->tmp_1x3_1;
  float* f_car = filter_->tmp_1x3_2;

  Transpose(R_wc, RwcT, 3, 3);
  Transpose(R_cm, RcmT, 3, 3);

  // Rwc^T * g
  MatMult(RwcT, kG, RwcTg, 3, 3, 1);

  float v_vec[3] = {0.f, velocity_cur, 0.f};
  float w_cross_v[3];
  // w x v:
  w_cross_v[0] = gyr_imu[1] * v_vec[2] - gyr_imu[2] * v_vec[1];
  w_cross_v[1] = gyr_imu[2] * v_vec[0] - gyr_imu[0] * v_vec[2];
  w_cross_v[2] = gyr_imu[0] * v_vec[1] - gyr_imu[1] * v_vec[0];

  // f_car = a_body - Rwc^T * g
  f_car[0] = w_cross_v[0] - RwcTg[0];
  f_car[1] = dv_dt + w_cross_v[1] - RwcTg[1];
  f_car[2] = w_cross_v[2] - RwcTg[2];

  MatMult(RcmT, f_car, h, 3, 3, 1);
}

void ComputeObservationJacobian(const float R_wc[9], const float R_cm[9],
                                float dv_dt, float velocity_cur,
                                const float gyr_imu[3], float H_3x6[18]) {
  float h_local[3];
  ComputeObservation(R_wc, R_cm, dv_dt, velocity_cur, gyr_imu, h_local);
  float* RwcT = filter_->tmp_3x3_3;
  float* RcmT = filter_->tmp_3x3_4;
  float* RwcTg = filter_->tmp_1x3_2;
  float* skew_tmp = filter_->tmp_3x3_5;
  float* H_left = filter_->tmp_3x3_1;
  float* tmp_3x3 = filter_->tmp_3x3_2;

  Transpose(R_wc, RwcT, 3, 3);
  Transpose(R_cm, RcmT, 3, 3);

  // d(h)/d(Rwc) = -Rcm^T * [Rwc^T * g]_x
  MatMult(RwcT, kG, RwcTg, 3, 3, 1);
  Skew3(RwcTg, skew_tmp);
  MatMult(RcmT, skew_tmp, H_left, 3, 3, 3);
  for (int i = 0; i < 9; ++i) H_left[i] = -H_left[i];

  // d(h)/d(Rcm) = [h]_x
  Skew3(h_local, tmp_3x3);

  for (int i = 0; i < 3; ++i) {
    const int i6 = i * 6, i3 = i * 3;
    std::memcpy(H_3x6 + i6, H_left + i3, 3 * sizeof(float));
    std::memcpy(H_3x6 + i6 + 3, tmp_3x3 + i3, 3 * sizeof(float));
  }
}

void InitFilter(const FilterParams& params) {
  g_params = params;
  if (filter_) std::free(filter_);
  filter_ = static_cast<filter*>(std::malloc(sizeof(filter)));
  MatSetIdentity(filter_->R_wc_, 3);
  MatSetIdentity(filter_->R_cm_, 3);
  std::memset(filter_->P_, 0, sizeof(filter_->P_));
  filter_->P_[0] = filter_->P_[7] = filter_->P_[14] = g_params.p0_rwc;
  filter_->P_[21] = g_params.p0_rcm_pitch;
  filter_->P_[28] = filter_->P_[35] = g_params.p0_rcm_roll_yaw;
  filter_->is_initialized_ = false;
  filter_->last_timestamp_ = 0;
  filter_->last_velocity_ = 0.f;
}

void Predict(float dt_s, const float gyr_imu[3]) {
  if (dt_s <= 0.f) return;
  float* theta = filter_->tmp_1x3_1;
  MatMult(filter_->R_cm_, gyr_imu, theta, 3, 3, 1);
  theta[0] *= dt_s;
  theta[1] *= dt_s;
  theta[2] *= dt_s;
  float* R_delta = filter_->tmp_3x3_1;
  RotVecToR(theta, R_delta);
  float* R_wc_new = filter_->tmp_3x3_2;
  MatMult(filter_->R_wc_, R_delta, R_wc_new, 3, 3, 3);
  std::memcpy(filter_->R_wc_, R_wc_new, 9 * sizeof(float));
  OrthonormalizeR3(filter_->R_wc_);
  const float qwc = g_params.q_wc, qcm = g_params.q_cm;
  float q2_dt = dt_s * qwc * qwc;
  filter_->P_[0] += q2_dt;
  filter_->P_[7] += q2_dt;
  filter_->P_[14] += q2_dt;
  filter_->P_[21] += dt_s * qcm * qcm;
  filter_->P_[28] += dt_s * qcm * qcm;
  filter_->P_[35] += dt_s * qcm * qcm;
}

void Update(int64_t timestamp, float velocity, const float acc_imu[3], const float gyr_imu[3]) {
  if (!filter_->is_initialized_) {
    // Initialization: just store velocity and timestamp
    filter_->last_timestamp_ = timestamp;
    filter_->last_velocity_ = velocity;
    filter_->is_initialized_ = true;
    return;
  }
  float dt_s = (timestamp - filter_->last_timestamp_) * 1e-9f;
  if (dt_s < 1e-6f) return;

  float dv_dt = (velocity - filter_->last_velocity_) / dt_s;

  const float max_dv = g_params.max_dv_dt_mps2;
  if (dv_dt > max_dv) dv_dt = max_dv;
  if (dv_dt < -max_dv) dv_dt = -max_dv;

  if (std::fabs(dv_dt) > g_params.dv_dt_skip_threshold) {
    filter_->last_timestamp_ = timestamp;
    filter_->last_velocity_ = velocity;
    return;
  }

  // Observation: h = Rcm^T * f_car
  float* h = filter_->tmp_1x3_3;
  ComputeObservation(filter_->R_wc_, filter_->R_cm_, dv_dt, velocity, gyr_imu, h);

  // Innovation: y = z - h(x)
  float* y = filter_->tmp_1x3_4;
  Vec3Sub(acc_imu, h, y);

  float* H = filter_->tmp_3x6_1;
  ComputeObservationJacobian(filter_->R_wc_, filter_->R_cm_, dv_dt, velocity, gyr_imu, H);

  // Ht -> tmp_3x6_2, PHt -> tmp_3x6_3
  float* Ht = filter_->tmp_3x6_2;
  Transpose(H, Ht, 3, 6);
  float* PHt = filter_->tmp_3x6_3;
  MatMult(filter_->P_, Ht, PHt, 6, 6, 3);

  float* S = filter_->tmp_3x3_3;
  MatMult(H, PHt, S, 3, 6, 3);

  float r2 = g_params.r_acc * g_params.r_acc;
  float* R_diag = filter_->tmp_3x3_4;
  std::memset(R_diag, 0, 9 * sizeof(float));
  R_diag[0] = R_diag[4] = R_diag[8] = r2;

  MatAdd(S, R_diag, filter_->tmp_3x3_5, 3, 3);
  float* Si = filter_->tmp_3x3_3;
  if (!Inv3x3(filter_->tmp_3x3_5, Si)) return;

  // K -> tmp_3x6_2
  float* K = filter_->tmp_3x6_2;
  MatMult(PHt, Si, K, 6, 3, 3);
  float* dx = filter_->tmp_1x6_1;
  MatMult(K, y, dx, 6, 3, 1);

  // State injection (right perturbation)
  float* R_delta = filter_->tmp_3x3_1;
  float* R_new = filter_->tmp_3x3_2;
  RotVecToR(dx, R_delta);
  MatMult(filter_->R_wc_, R_delta, R_new, 3, 3, 3);
  std::memcpy(filter_->R_wc_, R_new, 9 * sizeof(float));
  OrthonormalizeR3(filter_->R_wc_);
  RotVecToR(dx + 3, R_delta);
  MatMult(filter_->R_cm_, R_delta, R_new, 3, 3, 3);
  std::memcpy(filter_->R_cm_, R_new, 9 * sizeof(float));
  OrthonormalizeR3(filter_->R_cm_);

  // P = (I-KH)*P
  float* KH = filter_->tmp_6x6_1;
  MatMult(K, H, KH, 6, 3, 6);
  float* IKH = filter_->tmp_6x6_2;
  for (int i = 0; i < 36; ++i) IKH[i] = -KH[i];
  for (int i = 0; i < 6; ++i) IKH[i * 6 + i] += 1.f;
  float* P_new = filter_->tmp_6x6_3;
  MatMult(IKH, filter_->P_, P_new, 6, 6, 6);
  std::memcpy(filter_->P_, P_new, 36 * sizeof(float));

  // Store state for next step
  filter_->last_timestamp_ = timestamp;
  filter_->last_velocity_ = velocity;
}

void GetRcm(float* result) {
  std::memcpy(result, filter_->R_cm_, 9 * sizeof(float));
}

void GetRwc(float* result) {
  std::memcpy(result, filter_->R_wc_, 9 * sizeof(float));
}

}
