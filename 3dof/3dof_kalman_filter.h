#pragma once

#include <math.h>
#include <cstdint>
#include <cstring>

namespace mobili::ins {

// Nominal state: R_wc, R_cm (row-major 3x3 rotation matrices)
// Error state: 6-dim [δθ_wc(3), δθ_cm(3)] zero-mean, covariance P
// Predict: R_cm unchanged; R_wc updated by gyro integration; P plus process noise.
// Update: observation (temp.md model)
//   Using kinematic acceleration in body frame:
//   a_body = [0, (v_t - v_{t-1})/dt, 0] + w_gyro x [0, v_t, 0]
//   f_car = a_body - Rwc^T * g
//   h = Rcm^T * f_car
//
// ---------- Jacobian H (right perturbation) ----------
// Let: g = world-frame gravity [0,0,9.81]^T;
//      f_car = (v_t^C - Rwc^T*v_{t-1}^W)/dt - Rwc^T*g (X-right Y-forward Z-up);
//      h = Rcm^T * f_car.
//      Error state x = [δθ_wc; δθ_cm]
//      Right perturbation R' = R*(I+[δθ]_×) => (R')^T = (I-[δθ]_×)*R^T.
//
// (1) ∂h/∂δθ_cm:
//     h = Rcm^T*f_car;
//     Perturb Rcm^T => h' = (I-[δθ_cm]_×)*Rcm^T*f_car = h + [h]_×*δθ_cm;
//     thus ∂h/∂δθ_cm = [h]_×;
//
// (2) ∂h/∂δθ_wc:
//     f_car depends on Rwc^T*v_{t-1}^W and Rwc^T*g;
//     thus ∂f_car/∂δθ_wc = -[Rwc^T*g]_× - (1/dt)*[Rwc^T*v_{t-1}^W]_×,
//     hence ∂h/∂δθ_wc = Rcm^T * ∂f_car/∂δθ_wc
//                      = -Rcm^T*[Rwc^T*g]_× - (1/dt)*Rcm^T*[Rwc^T*v_{t-1}^W]_×.
//
// H = [ ∂h/∂δθ_wc | ∂h/∂δθ_cm ]
//   = [ -Rcm^T*[Rwc^T*g]_× - (1/dt)*Rcm^T*[Rwc^T*v_{t-1}^W]_× | [h]_× ]
struct filter {
  float R_wc_[9];  // vehicle-to-world rotation, row-major 3x3
  float R_cm_[9];  // device-to-vehicle rotation, row-major 3x3
  float P_[36];    // error state covariance 6x6, [δθ_wc(3), δθ_cm(3)]
  bool is_initialized_;
  int64_t last_timestamp_;
  float last_velocity_;

  // Temporary buffers for intermediate results, to avoid repeated allocation
  // for 3×3 matrices
  float tmp_3x3_1[9];
  float tmp_3x3_2[9];
  float tmp_3x3_3[9];
  float tmp_3x3_4[9];
  float tmp_3x3_5[9];
  // for 1×3 vectors
  float tmp_1x3_1[3];
  float tmp_1x3_2[3];
  float tmp_1x3_3[3];
  float tmp_1x3_4[3];
  // for 3×6 matrices
  float tmp_3x6_1[18];
  float tmp_3x6_2[18];
  float tmp_3x6_3[18];
  // for 6×6 matrices
  float tmp_6x6_1[36];
  float tmp_6x6_2[36];
  float tmp_6x6_3[36];
  // for 1×6 vectors
  float tmp_1x6_1[6];
  float tmp_1x6_2[6];
};

struct FilterParams {
  float q_wc = 0.01f;              // R_wc process noise [rad/√s]
  float q_cm = 1e-8f;              // R_cm process noise (roll/pitch/yaw)
  float r_acc = 0.08f;             // accelerometer observation noise [m/s²]
  float max_dv_dt_mps2 = 8.f;      // clamp |dv/dt| to boundary
  float p0_rwc = 1.f;
  float p0_rcm_pitch = 1e-8f;
  float p0_rcm_roll_yaw = 1e-8f;
  float dv_dt_skip_threshold = 1.0f;  // skip update when |dv/dt| > this [m/s²]
};

void InitFilter(const FilterParams& params = FilterParams());

void Predict(float dt_s, const float gyr_imu[3]);
void Update(int64_t timestamp, float velocity, const float acc_imu[3], const float gyr_imu[3]);
void GetRcm(float* result);
void GetRwc(float* result);
void SetDebugOutput(bool enable);

// Observation and Jacobian (dv_dt is precomputed by Update)
void ComputeObservation(const float R_wc[9], const float R_cm[9],
                        float dv_dt, float velocity_cur,
                        const float gyr_imu[3], float h[3]);
void ComputeObservationJacobian(const float R_wc[9], const float R_cm[9],
                                float dv_dt, float velocity_cur,
                                const float gyr_imu[3], float H_3x6[18]);

}
