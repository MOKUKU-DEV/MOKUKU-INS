#include "3dof_kalman_filter.h"
#include "3dof_kalman_filter_mat.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static float RotationAngleErrorRad(const float R_est[9], const float R_true[9]) {
  float R_err[9];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      R_err[i * 3 + j] = 0.f;
      for (int k = 0; k < 3; ++k)
        R_err[i * 3 + j] += R_est[k * 3 + i] * R_true[k * 3 + j];
    }
  }
  float trace = R_err[0] + R_err[4] + R_err[8];
  float cos_angle = (trace - 1.f) * 0.5f;
  if (cos_angle > 1.f) cos_angle = 1.f;
  if (cos_angle < -1.f) cos_angle = -1.f;
  return std::acos(cos_angle);
}

static const double kRadToDeg = 57.29577951308232;

struct Config {
  std::string data_dir = "/home/dmgz/ZWH/MOKUKU-INS/3dof/test_data";
  double start_ratio = 0.0;
  double end_ratio = 1.0;
  int update_interval = 1;
  bool inverse_imu_acc = true;
  double estimate_bias_seconds = 6.0;
  double gyr_bias_max_rad_s = 0.02;
  int vel_smooth_window = 5;
  double vel_motion_threshold = 0.1;

  // Filter parameters
  float q_wc = 0.01f;
  float q_cm = 1e-8f;
  float r_acc = 0.08f;
  float max_dv_dt_mps2 = 8.f;
  float p0_rwc = 1.f;
  float p0_rcm_pitch = 1e-8f;
  float p0_rcm_roll_yaw = 1e-8f;
  float dv_dt_skip_threshold = 1.0f;
};

static void Log(const char* level, const std::string& msg) {
  std::cerr << "[" << level << "] " << msg << std::endl;
}
static void LogInfo(const std::string& msg) { Log("INFO", msg); }
static void LogError(const std::string& msg) { Log("ERROR", msg); }

template <typename... Args>
static std::string StrCat(Args&&... args) {
  std::ostringstream oss;
  (oss << ... << std::forward<Args>(args));
  return oss.str();
}

static const float kGravity = 9.81f;

static float NormalizeAngleDiff(float rad) {
  const float kPi = 3.14159265358979323846f;
  rad = std::fmod(rad + kPi, 2.f * kPi);
  if (rad < 0.f) rad += 2.f * kPi;
  return rad - kPi;
}

static void ExtractYawPitchFromRwc(const float Rwc[9], float* yaw_rad, float* pitch_rad) {
  *yaw_rad = std::atan2(-Rwc[1], Rwc[4]);
  float r07 = Rwc[7];
  if (r07 > 1.f) r07 = 1.f;
  if (r07 < -1.f) r07 = -1.f;
  *pitch_rad = std::asin(-r07);
}

// Quaternion (qx,qy,qz,qw) -> row-major 3x3 rotation matrix R
static void QuatToRwc(float qx, float qy, float qz, float qw, float R[9]) {
  float qx2 = qx * qx, qy2 = qy * qy, qz2 = qz * qz;
  R[0] = 1.f - 2.f * (qy2 + qz2);
  R[1] = 2.f * (qx * qy - qz * qw);
  R[2] = 2.f * (qx * qz + qy * qw);
  R[3] = 2.f * (qx * qy + qz * qw);
  R[4] = 1.f - 2.f * (qx2 + qz2);
  R[5] = 2.f * (qy * qz - qx * qw);
  R[6] = 2.f * (qx * qz - qy * qw);
  R[7] = 2.f * (qy * qz + qx * qw);
  R[8] = 1.f - 2.f * (qx2 + qy2);
}

struct ImuSample {
  int64_t t_ns;
  float acc[3];
  float gyr[3];
};

struct TrajPose {
  int64_t t_ns;
  float Rwc[9];
  float vel;
};

// Parse IMU CSV: header "t_ns,acc_x,acc_y,acc_z,gyr_x,gyr_y,gyr_z"
static bool LoadImuCsv(const std::string& path, std::vector<ImuSample>* out) {
  std::ifstream file(path);
  if (!file.is_open()) {
    LogError("Failed to open IMU CSV: " + path);
    return false;
  }
  std::string header;
  if (!std::getline(file, header)) {
    LogError("Empty IMU file: " + path);
    return false;
  }
  out->clear();
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string cell;
    std::vector<double> row;
    while (std::getline(ss, cell, ',')) {
      try {
        row.push_back(std::stod(cell));
      } catch (...) {
        LogError("Invalid number in IMU line: " + cell);
        return false;
      }
    }
    if (row.size() != 7) {
      LogError("IMU line must have 7 columns, got " + std::to_string(row.size()));
      return false;
    }
    ImuSample s;
    s.t_ns = static_cast<int64_t>(row[0]);
    s.acc[0] = static_cast<float>(row[1]);
    s.acc[1] = static_cast<float>(row[2]);
    s.acc[2] = static_cast<float>(row[3]);
    s.gyr[0] = static_cast<float>(row[4]);
    s.gyr[1] = static_cast<float>(row[5]);
    s.gyr[2] = static_cast<float>(row[6]);
    out->push_back(s);
  }
  LogInfo("Loaded " + std::to_string(out->size()) + " IMU samples from " + path);
  return true;
}

// Parse trajectory: first line "origin_x origin_y origin_z", then "timestamp,px,py,pz,qx,qy,qz,qw,vel"
static bool LoadTrajectoryTxt(const std::string& path, float origin[3],
                              std::vector<TrajPose>* out) {
  std::ifstream file(path);
  if (!file.is_open()) {
    LogError(StrCat("Failed to open trajectory: ", path));
    return false;
  }
  std::string line0;
  if (!std::getline(file, line0)) {
    LogError("Empty trajectory file: " + path);
    return false;
  }
  std::stringstream ss0(line0);
  for (int i = 0; i < 3; ++i) {
    if (!(ss0 >> origin[i])) {
      LogError("Trajectory first line must have 3 doubles (origin)");
      return false;
    }
  }
  out->clear();
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string cell;
    std::vector<double> row;
    while (std::getline(ss, cell, ',')) {
      try {
        row.push_back(std::stod(cell));
      } catch (...) {
        LogError("Invalid number in trajectory line: " + cell);
        return false;
      }
    }
    if (row.size() != 9) {
      LogError(StrCat("Trajectory line must have 9 columns (t,px,py,pz,qx,qy,qz,qw,vel), got ", row.size()));
      return false;
    }
    TrajPose p;
    p.t_ns = static_cast<int64_t>(row[0]);
    QuatToRwc(static_cast<float>(row[4]), static_cast<float>(row[5]),
              static_cast<float>(row[6]), static_cast<float>(row[7]), p.Rwc);
    p.vel = static_cast<float>(row[8]);
    out->push_back(p);
  }
  LogInfo(StrCat("Loaded ", out->size(), " trajectory poses from ", path));
  return true;
}

// Convert trajectory rotations from absolute (ENU) to relative to first frame:
// R_relative_k = R_0^T * R_k (first frame becomes identity)
static void TrajectoryToRelativeFirstFrame(std::vector<TrajPose>* traj) {
  if (traj->empty() || traj->size() == 1) return;
  float R0_inv[9];
  mobili::ins::Transpose((*traj)[0].Rwc, R0_inv, 3, 3);
  float tmp[9];
  for (size_t k = 0; k < traj->size(); ++k) {
    mobili::ins::MatMult(R0_inv, (*traj)[k].Rwc, tmp, 3, 3, 3);
    for (int i = 0; i < 9; ++i) (*traj)[k].Rwc[i] = tmp[i];
  }
  LogInfo("Trajectory rotations converted to relative to first frame (R0^T * R_k).");
}

// Smooth trajectory scalar velocity with a simple moving average over samples (time-ordered).
// This reduces dv/dt spikes introduced by zero-velocity constraints or noisy velocity estimates.
static void SmoothTrajectoryVelocity(std::vector<TrajPose>* traj, int window_size) {
  if (traj->empty() || window_size <= 1) return;
  if (window_size % 2 == 0) ++window_size;  // ensure odd for symmetric window
  const int half = window_size / 2;
  const size_t N = traj->size();
  std::vector<float> orig_vel(N);
  for (size_t i = 0; i < N; ++i) orig_vel[i] = (*traj)[i].vel;
  for (size_t i = 0; i < N; ++i) {
    int i0 = static_cast<int>(i) - half;
    int i1 = static_cast<int>(i) + half;
    if (i0 < 0) i0 = 0;
    if (i1 >= static_cast<int>(N)) i1 = static_cast<int>(N) - 1;
    double sum = 0.0;
    int cnt = 0;
    for (int j = i0; j <= i1; ++j) {
      sum += static_cast<double>(orig_vel[static_cast<size_t>(j)]);
      ++cnt;
    }
    (*traj)[i].vel = static_cast<float>(sum / static_cast<double>(cnt));
  }
  LogInfo(StrCat("Trajectory velocity smoothed with window=", window_size, " samples (moving average)."));
}

// Interpolate velocity (linear in time); nearest pose for Rwc
static void InterpolateTrajectory(int64_t t_ns, const std::vector<TrajPose>& traj,
                                  float Rwc_out[9], float* vel_out) {
  if (traj.empty()) {
    mobili::ins::MatSetIdentity(Rwc_out, 3);
    *vel_out = 0.f;
    return;
  }
  if (traj.size() == 1) {
    for (int i = 0; i < 9; ++i) Rwc_out[i] = traj[0].Rwc[i];
    *vel_out = traj[0].vel;
    return;
  }
  // Nearest pose by time for Rwc
  size_t best = 0;
  int64_t best_dt = std::abs(static_cast<int64_t>(traj[0].t_ns - t_ns));
  for (size_t i = 1; i < traj.size(); ++i) {
    int64_t dt = std::abs(static_cast<int64_t>(traj[i].t_ns - t_ns));
    if (dt < best_dt) {
      best_dt = dt;
      best = i;
    }
  }
  for (int i = 0; i < 9; ++i) Rwc_out[i] = traj[best].Rwc[i];
  // Linear interpolation for velocity
  if (t_ns <= traj.front().t_ns) {
    *vel_out = traj.front().vel;
    return;
  }
  if (t_ns >= traj.back().t_ns) {
    *vel_out = traj.back().vel;
    return;
  }
  size_t i1 = 0;
  for (size_t i = 0; i + 1 < traj.size(); ++i) {
    if (traj[i].t_ns <= t_ns && t_ns <= traj[i + 1].t_ns) {
      i1 = i;
      break;
    }
  }
  int64_t t0 = traj[i1].t_ns, t1 = traj[i1 + 1].t_ns;
  double alpha = (t1 == t0) ? 1.0 : static_cast<double>(t_ns - t0) / static_cast<double>(t1 - t0);
  *vel_out = static_cast<float>((1.0 - alpha) * traj[i1].vel + alpha * traj[i1 + 1].vel);
}

int main() {
  Config cfg;
  std::string data_dir = cfg.data_dir;
  if (data_dir.back() != '/' && data_dir.back() != '\\') data_dir += '/';
  const std::string imu_csv_path = data_dir + "imu.csv";
  const std::string trajectory_txt_path = data_dir + "odom.txt";

  LogInfo(StrCat("Data dir: ", data_dir, " (imu.csv, odom.txt)"));

  std::vector<ImuSample> imu_samples;
  if (!LoadImuCsv(imu_csv_path, &imu_samples)) {
    return 1;
  }
  if (imu_samples.empty()) {
    LogError("No IMU samples.");
    return 1;
  }

  float origin[3];
  std::vector<TrajPose> traj;
  if (!LoadTrajectoryTxt(trajectory_txt_path, origin, &traj)) {
    return 1;
  }
  if (traj.empty()) {
    LogError("No trajectory poses.");
    return 1;
  }
  // Smooth velocity in ground-truth trajectory to reduce sharp steps from zero-velocity constraints, avoiding non-physical dv/dt spikes in observation model
  if (cfg.vel_smooth_window > 1) {
    SmoothTrajectoryVelocity(&traj, cfg.vel_smooth_window);
  }
  // Ground-truth trajectory is ENU absolute rotation; convert to relative rotation w.r.t. first frame
  TrajectoryToRelativeFirstFrame(&traj);

  // Use only IMU within trajectory time range [traj_min_ns, traj_max_ns]
  const int64_t traj_min_ns = traj.front().t_ns;
  const int64_t traj_max_ns = traj.back().t_ns;
  size_t imu_traj_start = 0;
  for (; imu_traj_start < imu_samples.size(); ++imu_traj_start) {
    if (imu_samples[imu_traj_start].t_ns >= traj_min_ns) break;
  }
  size_t imu_traj_end = imu_traj_start;
  for (size_t k = imu_traj_start; k < imu_samples.size(); ++k) {
    if (imu_samples[k].t_ns <= traj_max_ns) imu_traj_end = k;
    else break;
  }
  const size_t N_traj = (imu_traj_end >= imu_traj_start) ? (imu_traj_end - imu_traj_start + 1) : 0;
  LogInfo(StrCat("[Trajectory clip] Traj time ", (traj_min_ns * 1e-9), " ~ ", (traj_max_ns * 1e-9),
                " s -> IMU indices [", imu_traj_start, ", ", (imu_traj_end + 1), "), N=", N_traj));

  size_t start_idx = imu_traj_start + static_cast<size_t>(N_traj * cfg.start_ratio);
  size_t end_idx = imu_traj_start + static_cast<size_t>(N_traj * cfg.end_ratio);
  if (end_idx > imu_traj_end + 1) end_idx = imu_traj_end + 1;
  LogInfo(StrCat("[Ratio segment] start_ratio=", cfg.start_ratio, " end_ratio=", cfg.end_ratio,
                " over N_traj=", N_traj, " -> IMU indices [", start_idx, ", ", end_idx, ")"));
  if (start_idx >= end_idx) {
    LogError("Empty range: start_ratio/end_ratio (within trajectory clip).");
    return 1;
  }

  // Extrinsic/axis convention check: print first-frame IMU and pose for axis/sign alignment verification
  {
    const auto& im0 = imu_samples[start_idx];
    float acc_norm = std::sqrt(im0.acc[0] * im0.acc[0] + im0.acc[1] * im0.acc[1] + im0.acc[2] * im0.acc[2]);
    LogInfo(StrCat("[Convention check] First IMU (t=", (im0.t_ns * 1e-9), "s): acc=[",
                  im0.acc[0], ",", im0.acc[1], ",", im0.acc[2], "] |acc|=", acc_norm,
                  " (expect ~9.81 when level; up-axis should be ~-9.81, else try --inverse_imu_acc)"));
    LogInfo(StrCat("  gyr=[", im0.gyr[0], ",", im0.gyr[1], ",", im0.gyr[2], "]"));
    const auto& t0 = traj[0];
    float body_z[3] = {t0.Rwc[6], t0.Rwc[7], t0.Rwc[8]};  // body Z in world (row 2)
    LogInfo(StrCat("[Convention check] First pose (relative): body_Z_in_world=[", body_z[0], ",",
                  body_z[1], ",", body_z[2], "] (expect [0,0,1] when level), vel=", t0.vel));
  }

  const int update_interval = cfg.update_interval;  // 0 = no Update (predict only)

  // Bias: use first N seconds from segment start as bias reference
  float bias_acc[3] = {0.f, 0.f, 0.f};
  float bias_gyr[3] = {0.f, 0.f, 0.f};
  if (cfg.estimate_bias_seconds > 0 && start_idx < imu_samples.size()) {
    const int64_t t0_ns = imu_samples[start_idx].t_ns;
    const int64_t period_ns = static_cast<int64_t>(cfg.estimate_bias_seconds * 1e9);
    const float gyr_thresh = static_cast<float>(cfg.gyr_bias_max_rad_s);
    double sum_acc[3] = {0, 0, 0}, sum_gyr[3] = {0, 0, 0};
    int count = 0, count_gyr = 0;
    const float g_world[3] = {0.f, 0.f, kGravity};
    double sum_g_car[3] = {0, 0, 0};
    for (size_t i = start_idx; i < imu_samples.size() && (imu_samples[i].t_ns - t0_ns) <= period_ns; ++i) {
      for (int k = 0; k < 3; ++k) sum_acc[k] += imu_samples[i].acc[k];
      float gyr_norm =
          std::sqrt(imu_samples[i].gyr[0] * imu_samples[i].gyr[0] + imu_samples[i].gyr[1] * imu_samples[i].gyr[1] +
                    imu_samples[i].gyr[2] * imu_samples[i].gyr[2]);
      if (gyr_norm < gyr_thresh) {
        for (int k = 0; k < 3; ++k) sum_gyr[k] += imu_samples[i].gyr[k];
        ++count_gyr;
      }
      float Rwc[9], vel;
      InterpolateTrajectory(imu_samples[i].t_ns, traj, Rwc, &vel);
      float RwcT[9], g_car[3];
      mobili::ins::Transpose(Rwc, RwcT, 3, 3);
      mobili::ins::MatMult(RwcT, g_world, g_car, 3, 3, 1);
      for (int k = 0; k < 3; ++k) sum_g_car[k] += g_car[k];
      ++count;
    }
    if (count > 0) {
      const float sign = cfg.inverse_imu_acc ? -1.f : 1.f;
      for (int k = 0; k < 3; ++k) {
        bias_acc[k] = static_cast<float>(sum_acc[k] / count + sign * (sum_g_car[k] / count));
        bias_gyr[k] = (count_gyr > 0) ? static_cast<float>(sum_gyr[k] / count_gyr) : 0.f;
      }
      LogInfo(StrCat("[Bias estimate] ", cfg.estimate_bias_seconds, " s from segment start (t0=", (t0_ns * 1e-9),
                    "s): acc from ", count, " samples, gyr from ", count_gyr, " with |gyr|<",
                    cfg.gyr_bias_max_rad_s, " rad/s; bias_acc=[", bias_acc[0], ",", bias_acc[1], ",",
                    bias_acc[2], "] m/s^2, bias_gyr=[", bias_gyr[0], ",", bias_gyr[1], ",", bias_gyr[2], "] rad/s"));
    }
  }

  double sum_rwc_error_rad = 0.0;
  double max_rwc_error_rad = 0.0;
  double sum_pitch_error_rad = 0.0;
  double sum_total_angle_error_rad = 0.0;
  double sum_total_pitch_error_rad = 0.0;
  double sum_rcm_pitch_deg = 0.0;
  int num_compare = 0;
  const int64_t t0_ns = imu_samples[0].t_ns;  // Reference time origin (aligned with bias)
  bool motion_start_logged = false;
  float vel_prev_log = 0.f;
  int64_t t_prev_log_ns = imu_samples[start_idx].t_ns;  // Previous output frame time for dv_dt logging

  mobili::ins::FilterParams fp;
  fp.q_wc = cfg.q_wc;
  fp.q_cm = cfg.q_cm;
  fp.r_acc = cfg.r_acc;
  fp.max_dv_dt_mps2 = cfg.max_dv_dt_mps2;
  fp.p0_rwc = cfg.p0_rwc;
  fp.p0_rcm_pitch = cfg.p0_rcm_pitch;
  fp.p0_rcm_roll_yaw = cfg.p0_rcm_roll_yaw;
  fp.dv_dt_skip_threshold = cfg.dv_dt_skip_threshold;
  mobili::ins::InitFilter(fp);

  for (size_t i = start_idx; i < end_idx; ++i) {
    const ImuSample& imu = imu_samples[i];
    float Rwc_gt[9];
    float vel_gt;
    InterpolateTrajectory(imu.t_ns, traj, Rwc_gt, &vel_gt);
    double t_s = (imu.t_ns - t0_ns) * 1e-9;

    if (i > start_idx) {
      float dt_s = static_cast<float>(imu_samples[i].t_ns - imu_samples[i - 1].t_ns) * 1e-9f;
      float gyr_prev[3] = {
          imu_samples[i - 1].gyr[0] - bias_gyr[0],
          imu_samples[i - 1].gyr[1] - bias_gyr[1],
          imu_samples[i - 1].gyr[2] - bias_gyr[2],
      };
      mobili::ins::Predict(dt_s, gyr_prev);
    }

    // Every update_interval frames: Update (interval=0: no Update)
    if (update_interval > 0 && static_cast<int>(i - start_idx + 1) % update_interval == 0) {
      float acc_avg[3], gyr_avg[3];
      if (update_interval == 1) {
        acc_avg[0] = imu.acc[0] - bias_acc[0];
        acc_avg[1] = imu.acc[1] - bias_acc[1];
        acc_avg[2] = imu.acc[2] - bias_acc[2];
        if (cfg.inverse_imu_acc) {
          acc_avg[0] = -acc_avg[0];
          acc_avg[1] = -acc_avg[1];
          acc_avg[2] = -acc_avg[2];
        }
        gyr_avg[0] = imu.gyr[0] - bias_gyr[0];
        gyr_avg[1] = imu.gyr[1] - bias_gyr[1];
        gyr_avg[2] = imu.gyr[2] - bias_gyr[2];
      } else {
        const size_t interval_start = i + 1 - static_cast<size_t>(update_interval);
        const float T = static_cast<float>(imu_samples[i].t_ns - imu_samples[interval_start].t_ns) * 1e-9f;
        if (T < 1e-9f) {
          acc_avg[0] = imu.acc[0] - bias_acc[0];
          acc_avg[1] = imu.acc[1] - bias_acc[1];
          acc_avg[2] = imu.acc[2] - bias_acc[2];
          if (cfg.inverse_imu_acc) {
            acc_avg[0] = -acc_avg[0];
            acc_avg[1] = -acc_avg[1];
            acc_avg[2] = -acc_avg[2];
          }
          gyr_avg[0] = imu.gyr[0] - bias_gyr[0];
          gyr_avg[1] = imu.gyr[1] - bias_gyr[1];
          gyr_avg[2] = imu.gyr[2] - bias_gyr[2];
        } else {
          // Preintegration: ΔR=I, Δv=0; recurse ΔR_{k+1}=ΔR_k·exp(ω_k dt), Δv_{k+1}=Δv_k+ΔR_k·a_k·dt; ω_eq=log(ΔR)/T, a_eq_end=ΔR^T·(Δv/T)
          float delta_R[9], delta_v[3];
          mobili::ins::MatSetIdentity(delta_R, 3);
          delta_v[0] = delta_v[1] = delta_v[2] = 0.f;
          float theta_k[3], exp_k[9], delta_R_a[3];
          for (size_t k = interval_start; k < i; ++k) {
            const float dt_k = static_cast<float>(imu_samples[k + 1].t_ns - imu_samples[k].t_ns) * 1e-9f;
            float gyr_k[3] = {
                imu_samples[k].gyr[0] - bias_gyr[0],
                imu_samples[k].gyr[1] - bias_gyr[1],
                imu_samples[k].gyr[2] - bias_gyr[2],
            };
            float acc_k[3] = {
                imu_samples[k].acc[0] - bias_acc[0],
                imu_samples[k].acc[1] - bias_acc[1],
                imu_samples[k].acc[2] - bias_acc[2],
            };
            if (cfg.inverse_imu_acc) {
              acc_k[0] = -acc_k[0];
              acc_k[1] = -acc_k[1];
              acc_k[2] = -acc_k[2];
            }
            for (int d = 0; d < 3; ++d) theta_k[d] = gyr_k[d] * dt_k;
            mobili::ins::RotVecToR(theta_k, exp_k);
            acc_k[0] *= dt_k;
            acc_k[1] *= dt_k;
            acc_k[2] *= dt_k;
            mobili::ins::MatMult(delta_R, acc_k, delta_R_a, 3, 3, 1);
            for (int d = 0; d < 3; ++d) delta_v[d] += delta_R_a[d];
            float delta_R_new[9];
            mobili::ins::MatMult(delta_R, exp_k, delta_R_new, 3, 3, 3);
            for (int d = 0; d < 9; ++d) delta_R[d] = delta_R_new[d];
          }
          float log_delta_R[3];
          mobili::ins::RToRotVec(delta_R, log_delta_R);
          for (int d = 0; d < 3; ++d) gyr_avg[d] = log_delta_R[d] / T;
          float delta_RT[9], a_eq_start[3];
          mobili::ins::Transpose(delta_R, delta_RT, 3, 3);
          for (int d = 0; d < 3; ++d) a_eq_start[d] = delta_v[d] / T;
          mobili::ins::MatMult(delta_RT, a_eq_start, acc_avg, 3, 3, 1);
        }
      }
      mobili::ins::Update(imu.t_ns, vel_gt, acc_avg, gyr_avg);
    }

    if (i % 2 == 0) {
      float Rwc_est[9], Rcm_est[9];
      mobili::ins::GetRwc(Rwc_est);
      mobili::ins::GetRcm(Rcm_est);
      float rwc_err_rad = RotationAngleErrorRad(Rwc_est, Rwc_gt);
      sum_rwc_error_rad += static_cast<double>(rwc_err_rad);
      if (rwc_err_rad > max_rwc_error_rad) max_rwc_error_rad = static_cast<double>(rwc_err_rad);
      float yaw_est, pitch_est, yaw_gt, pitch_gt;
      ExtractYawPitchFromRwc(Rwc_est, &yaw_est, &pitch_est);
      ExtractYawPitchFromRwc(Rwc_gt, &yaw_gt, &pitch_gt);
      float pitch_err_deg = NormalizeAngleDiff(pitch_est - pitch_gt) * kRadToDeg;
      sum_pitch_error_rad += static_cast<double>(std::fabs(NormalizeAngleDiff(pitch_est - pitch_gt)));

      // R_total = R_wc * R_cm (IMU-to-world); if IMU=car, should match Rwc_gt; diagnostic for R_cm absorbing pitch
      float R_total[9];
      mobili::ins::MatMult(Rwc_est, Rcm_est, R_total, 3, 3, 3);
      float total_err_rad = RotationAngleErrorRad(R_total, Rwc_gt);
      sum_total_angle_error_rad += static_cast<double>(total_err_rad);
      float yaw_tot, pitch_tot;
      ExtractYawPitchFromRwc(R_total, &yaw_tot, &pitch_tot);
      sum_total_pitch_error_rad += static_cast<double>(std::fabs(NormalizeAngleDiff(pitch_tot - pitch_gt)));
      float rcm_yaw, rcm_pitch;
      ExtractYawPitchFromRwc(Rcm_est, &rcm_yaw, &rcm_pitch);
      sum_rcm_pitch_deg += static_cast<double>(rcm_pitch * kRadToDeg);

      // Motion start: log once when velocity first exceeds threshold
      if (vel_gt >= static_cast<float>(cfg.vel_motion_threshold) && !motion_start_logged) {
        LogInfo(StrCat("[Motion start] frame=", i, ", t_s=", t_s, ", vel=", vel_gt,
                      " m/s, pitch_est=", (pitch_est * kRadToDeg), " deg, pitch_gt=", (pitch_gt * kRadToDeg),
                      " deg, pitch_err=", pitch_err_deg, " deg"));
        motion_start_logged = true;
      }

      // dv_dt (between output frames) for acceleration/deceleration and observation (from trajectory velocity diff)
      float dv_dt_log = 0.f;
      if (i > start_idx) {
        float dt_log_s = static_cast<float>(imu.t_ns - t_prev_log_ns) * 1e-9f;
        if (dt_log_s > 1e-9f) dv_dt_log = (vel_gt - vel_prev_log) / dt_log_s;
      }
      t_prev_log_ns = imu.t_ns;
      vel_prev_log = vel_gt;

      // Compute longitudinal acceleration a_y_imu from IMU specific force + gravity, compare with dv_dt:
      // In Update, acc_imu = f_car = a_body - Rwc^T*g, and a_body_y ≈ dv_dt (w×v y-component is 0 in current model).
      // So a_y_imu ≈ acc_imu_y + (Rwc^T*g)_y, should match dv_dt numerically.
      float acc_corr[3] = {
          imu.acc[0] - bias_acc[0],
          imu.acc[1] - bias_acc[1],
          imu.acc[2] - bias_acc[2],
      };
      if (cfg.inverse_imu_acc) {
        acc_corr[0] = -acc_corr[0];
        acc_corr[1] = -acc_corr[1];
        acc_corr[2] = -acc_corr[2];
      }
      const float g_world[3] = {0.f, 0.f, kGravity};
      float RwcT_gt[9], g_car_gt[3];
      mobili::ins::Transpose(Rwc_gt, RwcT_gt, 3, 3);
      mobili::ins::MatMult(RwcT_gt, g_world, g_car_gt, 3, 3, 1);
      float a_y_imu = acc_corr[1] + g_car_gt[1];

      // Pitch from raw IMU acc Y/Z (no gravity subtraction, sign correction only):
      // When level, acc may be [0,0,+g] or [0,0,-g] (inverse_imu_acc + sensor convention).
      // Normalize to: level 0°, nose-up positive.
      float acc_yz[3] = {
          imu.acc[0],
          imu.acc[1],
          imu.acc[2],
      };
      if (cfg.inverse_imu_acc) {
        acc_yz[0] = -acc_yz[0];
        acc_yz[1] = -acc_yz[1];
        acc_yz[2] = -acc_yz[2];
      }
      // Convention-independent: level 0°, nose-up positive. pitch=atan2(-acc_Y,acc_Z) if acc_z>=0, else atan2(acc_Y,-acc_Z)
      float pitch_from_acc_rad =
          (acc_yz[2] >= 0.f) ? std::atan2(-acc_yz[1], acc_yz[2]) : std::atan2(acc_yz[1], -acc_yz[2]);
      float pitch_from_acc_deg = pitch_from_acc_rad * static_cast<float>(kRadToDeg);

      ++num_compare;

      LogInfo(StrCat("Frame ", i, " t_s=", t_s, " vel=", vel_gt, " m/s dv_dt=", dv_dt_log,
                    " m/s^2 a_y_imu=", a_y_imu, " m/s^2 | pitch_est=", (pitch_est * kRadToDeg),
                    " deg pitch_gt=", (pitch_gt * kRadToDeg), " deg pitch_from_acc=", pitch_from_acc_deg,
                    " deg pitch_err=", pitch_err_deg, " deg | Rwc_angle_err=", (rwc_err_rad * kRadToDeg), " deg"));
    }
  }

  if (num_compare > 0) {
    double mean_rwc_rad = sum_rwc_error_rad / num_compare;
    double mean_pitch_rad = sum_pitch_error_rad / num_compare;
    double mean_total_rad = sum_total_angle_error_rad / num_compare;
    double mean_total_pitch_rad = sum_total_pitch_error_rad / num_compare;
    double mean_rcm_pitch_deg = sum_rcm_pitch_deg / num_compare;
    LogInfo(StrCat("[Summary] Rwc: mean angle error = ", mean_rwc_rad, " rad (",
                  (mean_rwc_rad * kRadToDeg), " deg), max = ", max_rwc_error_rad, " rad (",
                  (max_rwc_error_rad * kRadToDeg), " deg), over ", num_compare, " frames."));
    LogInfo(StrCat("  Rwc pitch mean error = ", (mean_pitch_rad * kRadToDeg), " deg."));
    LogInfo(StrCat("  R_total (=Rwc*Rcm, IMU-to-world) vs gt: mean angle = ", (mean_total_rad * kRadToDeg),
                  " deg, mean pitch err = ", (mean_total_pitch_rad * kRadToDeg),
                  " deg (if IMU=car and total<<Rwc err, pitch is in R_cm)."));
    LogInfo(StrCat("  R_cm mean pitch = ", mean_rcm_pitch_deg, " deg."));
  }

  LogInfo("Finished.");
  return 0;
}
