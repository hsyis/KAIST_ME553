#pragma once

struct Anymal {
    struct Joint {
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
        Eigen::Vector3d axis;
    };

    Joint base;
    Joint hip;
    Joint thigh;
    Joint shank;
    Joint adapter; // shank_to_adapter
    Joint foot;    // adapter_to_foot

    Anymal() {
        base.rpy = {0.0, 0.0, 0.0};
        base.xyz = {0.0, 0.0, 0.0};
        base.axis = {0, 0, 0};

        hip.rpy = {0.0, 0.0, 0.0};
        hip.xyz = {0.277, 0.116, 0.0};
        hip.axis = {1, 0, 0};

        thigh.rpy = {0.0, 0.0, 0.0};
        thigh.xyz = {0.0635, 0.041, 0.0};
        thigh.axis = {0, 1, 0};

        shank.rpy = {0.0, 0.0, 0.0};
        shank.xyz = {0.0, 0.109, -0.25};
        shank.axis = {0, 1, 0};

        adapter.rpy = {0.0, 0.0, 0.0};
        adapter.xyz = {0.1, -0.02, 0.0};
        adapter.axis = {0, 0, 0};

        foot.rpy = {0.0, 0.0, 0.0};
        foot.xyz = {0.0, 0.0, -0.32125};
        foot.axis = {0, 0, 0};
    }
};

Eigen::Matrix3d toRot(double w, double x, double y, double z) {
    Eigen::Matrix3d R;

    auto op1 = [](auto a, auto b) {
        return 1 - 2 * (a * a + b * b);
    };

    auto op2 = [](auto a, auto b, auto c, auto d) {
        return 2 * (a * b + c * d);
    };

    R << op1(y, z), op2(x, y, z, -w), op2(x, z, y, w),
         op2(x, y, z, w), op1(x, z), op2(y, z, x, -w),
         op2(x, z, y, -w), op2(y, z, x, w), op1(x, y);

    return R;
}

Eigen::Matrix3d toRot(const Eigen::VectorXd& rpy) {
    Eigen::Matrix3d R_x;
    {
        double c0 = std::cos(rpy[0]);
        double s0 = std::sin(rpy[0]);

        R_x << 1, 0, 0,
               0, c0, -s0,
               0, s0, c0;
    }

    Eigen::Matrix3d R_y;
    {
        double c1 = std::cos(rpy[1]);
        double s1 = std::sin(rpy[1]);

        R_y << c1, 0, s1,
               0, 1, 0,
               -s1, 0, c1;
    }

    Eigen::Matrix3d R_z;
    {
        double c2 = std::cos(rpy[2]);
        double s2 = std::sin(rpy[2]);

        R_z << c2, -s2, 0,
               s2, c2, 0,
               0, 0, 1;
    }

    return R_z * R_y * R_x;
}

enum class Rotation { roll, pitch, yaw };

Eigen::Matrix3d toRot(double theta, Rotation rotation) {
    Eigen::VectorXd rpy(3);

    switch (rotation) {
    case Rotation::roll:
        rpy << theta, 0, 0;
        break;
    case Rotation::pitch:
        rpy << 0, theta, 0;
        break;
    case Rotation::yaw:
        rpy << 0, 0, theta;
        break;
    }

    return toRot(rpy);
}

Eigen::Matrix3d skewSymmetric(const Eigen::Vector3d& a) {
    Eigen::Matrix3d m;
    m << 0, -a[2], a[1],
         a[2], 0, -a[0],
         -a[1], a[0], 0;
    return m;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Anymal anymal;

    Eigen::Matrix3d R_w_b = toRot(anymal.base.rpy);
    Eigen::Matrix3d R_b_1 = toRot(anymal.hip.rpy);
    Eigen::Matrix3d R_1_2 = toRot(anymal.thigh.rpy);
    Eigen::Matrix3d R_2_3 = toRot(anymal.shank.rpy);
    Eigen::Matrix3d R_3_4 = toRot(anymal.adapter.rpy);

    Eigen::Matrix3d R_b_b = toRot(gc[3], gc[4], gc[5], gc[6]);
    Eigen::Matrix3d R_1_1 = toRot(gc[7], Rotation::roll);
    Eigen::Matrix3d R_2_2 = toRot(gc[8], Rotation::pitch);
    Eigen::Matrix3d R_3_3 = toRot(gc[9], Rotation::pitch);

    Eigen::Matrix3d R_w_1 = R_w_b * R_b_b * R_b_1 * R_1_1;
    Eigen::Matrix3d R_w_2 = R_w_1 * R_1_2 * R_2_2;
    Eigen::Matrix3d R_w_3 = R_w_2 * R_2_3 * R_3_3;
    Eigen::Matrix3d R_w_4 = R_w_3 * R_3_4;

    Eigen::Vector3d r_b_1 = R_w_b * R_b_b * anymal.hip.xyz;
    Eigen::Vector3d r_1_2 = R_w_1 * anymal.thigh.xyz;
    Eigen::Vector3d r_2_3 = R_w_2 * anymal.shank.xyz;
    Eigen::Vector3d r_3_4 = R_w_3 * anymal.adapter.xyz;
    Eigen::Vector3d r_4_e = R_w_4 * anymal.foot.xyz;

    Eigen::Vector3d r_3_e = r_3_4 + r_4_e;
    Eigen::Vector3d r_2_e = r_2_3 + r_3_e;
    Eigen::Vector3d r_1_e = r_1_2 + r_2_e;
    Eigen::Vector3d r_b_e = r_b_1 + r_1_e;

    Eigen::Vector3d P_w_1 = R_w_1 * anymal.hip.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * anymal.thigh.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * anymal.shank.axis;

    Eigen::MatrixXd positional_jacobian(3, 18);

    positional_jacobian << Eigen::Matrix3d::Identity(),
                           -skewSymmetric(r_b_e),
                           -skewSymmetric(r_1_e) * P_w_1,
                           -skewSymmetric(r_2_e) * P_w_2,
                           -skewSymmetric(r_3_e) * P_w_3,
                           Eigen::MatrixXd::Zero(3, 9); // only LF

    return positional_jacobian * gv;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Anymal anymal;

    Eigen::Matrix3d R_w_b = toRot(anymal.base.rpy);
    Eigen::Matrix3d R_b_1 = toRot(anymal.hip.rpy);
    Eigen::Matrix3d R_1_2 = toRot(anymal.thigh.rpy);
    Eigen::Matrix3d R_2_3 = toRot(anymal.shank.rpy);

    Eigen::Matrix3d R_b_b = toRot(gc[3], gc[4], gc[5], gc[6]);
    Eigen::Matrix3d R_1_1 = toRot(gc[7], Rotation::roll);
    Eigen::Matrix3d R_2_2 = toRot(gc[8], Rotation::pitch);
    Eigen::Matrix3d R_3_3 = toRot(gc[9], Rotation::pitch);

    Eigen::Matrix3d R_w_1 = R_w_b * R_b_b * R_b_1 * R_1_1;
    Eigen::Matrix3d R_w_2 = R_w_1 * R_1_2 * R_2_2;
    Eigen::Matrix3d R_w_3 = R_w_2 * R_2_3 * R_3_3;

    Eigen::Vector3d P_w_1 = R_w_1 * anymal.hip.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * anymal.thigh.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * anymal.shank.axis;

    Eigen::MatrixXd angular_jacobian(3, 18);

    angular_jacobian << Eigen::Matrix3d::Zero(),
                        Eigen::Matrix3d::Identity(),
                        P_w_1,
                        P_w_2,
                        P_w_3,
                        Eigen::MatrixXd::Zero(3, 9); // only LF

    return angular_jacobian * gv;
}
