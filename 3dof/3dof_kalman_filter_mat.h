#pragma once

#include <math.h>
#include <cstring>

namespace mobili::ins {

// All matrix operations use row-major layout.

// At = A^T
// A is rows x cols, At is cols x rows
inline void Transpose(const float* A, float* At, int rows, int cols) {
  for (int i = 0; i < rows; ++i)
    for (int j = 0; j < cols; ++j)
      At[j * rows + i] = A[i * cols + j];
}

// C = A + B, A, B, C are rows x cols
inline void MatAdd(const float* A, const float* B, float* C, int rows, int cols) {
  for (int i = 0; i < rows * cols; ++i)
    C[i] = A[i] + B[i];
}

// C = A * B, A is arows x acols, B is acols x bcols, C is arows x bcols
inline void MatMult(const float* A, const float* B, float* C,
                    int arows, int acols, int bcols) {
  for (int i = 0; i < arows; ++i) {
    for (int j = 0; j < bcols; ++j) {
      C[i * bcols + j] = 0.f;
      for (int k = 0; k < acols; ++k)
        C[i * bcols + j] += A[i * acols + k] * B[k * bcols + j];
    }
  }
}

inline void MatSetIdentity(float* A, int n) {
  std::memset(A, 0, static_cast<size_t>(n * n) * sizeof(float));
  for (int i = 0; i < n; ++i) A[i * n + i] = 1.f;
}

inline void Vec3Sub(const float a[3], const float b[3], float out[3]) {
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}

// Skew-symmetric matrix [v]_x, [v]_x * u = v x u (cross product)
inline void Skew3(const float v[3], float S[9]) {
  S[0] = 0.f;
  S[1] = -v[2];
  S[2] = v[1];
  S[3] = v[2];
  S[4] = 0.f;
  S[5] = -v[0];
  S[6] = -v[1];
  S[7] = v[0];
  S[8] = 0.f;
}

// theta (axis-angle, |theta|=angle) -> rotation matrix R
// Rodrigues: R = I + (sin θ/θ)[θ]_× + ((1-cos θ)/θ²)(θθ^T)
inline void RotVecToR(const float theta[3], float R[9]) {
  float angle_sq = theta[0] * theta[0] + theta[1] * theta[1] + theta[2] * theta[2];
  if (angle_sq < 1e-14f) {
    MatSetIdentity(R, 3);
    return;
  }
  float angle = sqrtf(angle_sq);
  float co = cosf(angle), si = sinf(angle);
  float scale1 = si / angle;
  float scale2 = (1.f - co) / (angle_sq);
  float t0 = theta[0], t1 = theta[1], t2 = theta[2];
  R[0] = co + scale2 * t0 * t0;
  R[1] = -scale1 * t2 + scale2 * t0 * t1;
  R[2] = scale1 * t1 + scale2 * t0 * t2;
  R[3] = scale1 * t2 + scale2 * t1 * t0;
  R[4] = co + scale2 * t1 * t1;
  R[5] = -scale1 * t0 + scale2 * t1 * t2;
  R[6] = -scale1 * t1 + scale2 * t2 * t0;
  R[7] = scale1 * t0 + scale2 * t2 * t1;
  R[8] = co + scale2 * t2 * t2;
}

// Log map SO(3) -> so(3): R (row-major) -> rotation vector theta (rad), R = exp([theta]_x)
inline void RToRotVec(const float R[9], float theta[3]) {
  float tr = R[0] + R[4] + R[8];
  float cos_angle = (tr - 1.f) * 0.5f;
  if (cos_angle > 1.f) cos_angle = 1.f;
  if (cos_angle < -1.f) cos_angle = -1.f;
  float angle = acosf(cos_angle);
  const float kSinEps = 1e-6f;
  float sin_angle = sinf(angle);
  if (angle < kSinEps) {
    theta[0] = (R[7] - R[5]) * 0.5f;
    theta[1] = (R[2] - R[6]) * 0.5f;
    theta[2] = (R[3] - R[1]) * 0.5f;
    return;
  }
  if (sin_angle < kSinEps) {
    float ax = sqrtf((R[0] + 1.f) * 0.5f > 0.f ? (R[0] + 1.f) * 0.5f : 0.f);
    float ay = sqrtf((R[4] + 1.f) * 0.5f > 0.f ? (R[4] + 1.f) * 0.5f : 0.f);
    float az = sqrtf((R[8] + 1.f) * 0.5f > 0.f ? (R[8] + 1.f) * 0.5f : 0.f);
    if (ax >= ay && ax >= az && ax > kSinEps) {
      theta[0] = angle * ax;
      theta[1] = angle * (R[3] / (2.f * ax));
      theta[2] = angle * (R[6] / (2.f * ax));
    } else if (ay > kSinEps) {
      theta[0] = angle * (R[1] / (2.f * ay));
      theta[1] = angle * ay;
      theta[2] = angle * (R[7] / (2.f * ay));
    } else {
      theta[0] = angle * (R[2] / (2.f * az));
      theta[1] = angle * (R[5] / (2.f * az));
      theta[2] = angle * az;
    }
    return;
  }
  float scale = angle / (2.f * sin_angle);
  theta[0] = scale * (R[7] - R[5]);
  theta[1] = scale * (R[2] - R[6]);
  theta[2] = scale * (R[3] - R[1]);
}

// 3x3 row-major matrix orthonormalized by column Gram–Schmidt, keeping R in SO(3)
// Col0 normalized, col2 = col0 × col1 normalized, col1 normalized by col2 × col0
inline void OrthonormalizeR3(float R[9]) {
  float* c0 = &R[0];   // col0: R[0], R[3], R[6]
  float* c1 = &R[1];   // col1: R[1], R[4], R[7]
  float* c2 = &R[2];   // col2: R[2], R[5], R[8]
  float n0 = sqrtf(c0[0] * c0[0] + c0[3] * c0[3] + c0[6] * c0[6]);
  if (n0 < 1e-12f) return;
  n0 = 1.f / n0;
  c0[0] *= n0; c0[3] *= n0; c0[6] *= n0;
  float cx = c0[3] * c1[6] - c0[6] * c1[3];
  float cy = c0[6] * c1[0] - c0[0] * c1[6];
  float cz = c0[0] * c1[3] - c0[3] * c1[0];
  float n2 = sqrtf(cx * cx + cy * cy + cz * cz);
  if (n2 < 1e-12f) return;
  n2 = 1.f / n2;
  c2[0] = cx * n2; c2[3] = cy * n2; c2[6] = cz * n2;
  float bx = c2[3] * c0[6] - c2[6] * c0[3];
  float by = c2[6] * c0[0] - c2[0] * c0[6];
  float bz = c2[0] * c0[3] - c2[3] * c0[0];
  float n1 = sqrtf(bx * bx + by * by + bz * bz);
  if (n1 < 1e-12f) return;
  n1 = 1.f / n1;
  c1[0] = bx * n1; c1[3] = by * n1; c1[6] = bz * n1;
}

// 3x3 matrix inverse; returns false if singular, true otherwise
inline bool Inv3x3(const float S[9], float Si[9]) {
  float S00 = S[0], S01 = S[1], S02 = S[2];
  float S10 = S[3], S11 = S[4], S12 = S[5];
  float S20 = S[6], S21 = S[7], S22 = S[8];
  float det = S00 * (S11 * S22 - S12 * S21) -
              S01 * (S10 * S22 - S12 * S20) +
              S02 * (S10 * S21 - S11 * S20);
  if (det * det < 1e-24f) return false;
  float inv_det = 1.f / det;
  Si[0] = (S11 * S22 - S12 * S21) * inv_det;
  Si[1] = (S02 * S21 - S01 * S22) * inv_det;
  Si[2] = (S01 * S12 - S02 * S11) * inv_det;
  Si[3] = (S12 * S20 - S10 * S22) * inv_det;
  Si[4] = (S00 * S22 - S02 * S20) * inv_det;
  Si[5] = (S02 * S10 - S00 * S12) * inv_det;
  Si[6] = (S10 * S21 - S11 * S20) * inv_det;
  Si[7] = (S01 * S20 - S00 * S21) * inv_det;
  Si[8] = (S00 * S11 - S01 * S10) * inv_det;
  return true;
}

}
