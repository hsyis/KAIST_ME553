#pragma once

struct Kinova {
    struct Joint {
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
        Eigen::Vector3d axis;
    };

    Joint base;
    Joint joint1;
    Joint joint2;
    Joint joint3;
    Joint joint4;
    Joint joint5;
    Joint joint6;
    Joint end;

    Kinova() {
        base.rpy = {0, 0, 0};
        base.xyz = {0, 0, 0};
        base.axis = {0, 0, 0};

        joint1.rpy = {0, 3.14159265359, 0};
        joint1.xyz = {0, 0, 0.15675};
        joint1.axis = {0, 0, 1};

        joint2.rpy = {-1.57079632679, 0, 3.14159265359};
        joint2.xyz = {0, 0.0016, -0.11875};
        joint2.axis = {0, 0, 1};

        joint3.rpy = {0, 3.14159265359, 0};
        joint3.xyz = {0, -0.410, 0};
        joint3.axis = {0, 0, 1};

        joint4.rpy = {-1.57079632679, 0, 3.14159265359};
        joint4.xyz = {0, 0.2073, -0.0114};
        joint4.axis = {0, 0, 1};

        joint5.rpy = {1.57079632679, 0, 3.14159265359};
        joint5.xyz = {0, 0, -0.10375};
        joint5.axis = {0, 0, 1};

        joint6.rpy = {-1.57079632679, 0, 3.14159265359};
        joint6.xyz = {0, 0.10375, 0};
        joint6.axis = {0, 0, 1};

        end.rpy = {3.14159265359, 0, 0};
        end.xyz = {0, 0, -0.1600};
        end.axis = {0, 0, 0};
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

Eigen::Vector3d toAngleAxis(const Eigen::Matrix3d& rotation) {
    Eigen::AngleAxisd aa(rotation);
    return aa.axis() * aa.angle();
}

/// do not change the name of the method
inline Eigen::VectorXd getVelocityCommand (const Eigen::VectorXd& gc, const Eigen::Vector3d& pos, const Eigen::Vector4d& quat) {
    Kinova kinova;

    Eigen::Matrix3d R_w_b = toRot(kinova.base.rpy);
    Eigen::Matrix3d R_b_1 = toRot(kinova.joint1.rpy);
    Eigen::Matrix3d R_1_2 = toRot(kinova.joint2.rpy);
    Eigen::Matrix3d R_2_3 = toRot(kinova.joint3.rpy);
    Eigen::Matrix3d R_3_4 = toRot(kinova.joint4.rpy);
    Eigen::Matrix3d R_4_5 = toRot(kinova.joint5.rpy);
    Eigen::Matrix3d R_5_6 = toRot(kinova.joint6.rpy);
    Eigen::Matrix3d R_6_e = toRot(kinova.end.rpy);

    Eigen::Matrix3d R_1_1 = toRot(gc[0], Rotation::yaw);
    Eigen::Matrix3d R_2_2 = toRot(gc[1], Rotation::yaw);
    Eigen::Matrix3d R_3_3 = toRot(gc[2], Rotation::yaw);
    Eigen::Matrix3d R_4_4 = toRot(gc[3], Rotation::yaw);
    Eigen::Matrix3d R_5_5 = toRot(gc[4], Rotation::yaw);
    Eigen::Matrix3d R_6_6 = toRot(gc[5], Rotation::yaw);

    Eigen::Matrix3d R_w_1 = R_w_b * R_b_1 * R_1_1;
    Eigen::Matrix3d R_w_2 = R_w_1 * R_1_2 * R_2_2;
    Eigen::Matrix3d R_w_3 = R_w_2 * R_2_3 * R_3_3;
    Eigen::Matrix3d R_w_4 = R_w_3 * R_3_4 * R_4_4;
    Eigen::Matrix3d R_w_5 = R_w_4 * R_4_5 * R_5_5;
    Eigen::Matrix3d R_w_6 = R_w_5 * R_5_6 * R_6_6;
    Eigen::Matrix3d R_w_e = R_w_6 * R_6_e;

    Eigen::Vector3d r_w_b = kinova.base.xyz;
    Eigen::Vector3d r_b_1 = R_w_b * kinova.joint1.xyz;
    Eigen::Vector3d r_1_2 = R_w_1 * kinova.joint2.xyz;
    Eigen::Vector3d r_2_3 = R_w_2 * kinova.joint3.xyz;
    Eigen::Vector3d r_3_4 = R_w_3 * kinova.joint4.xyz;
    Eigen::Vector3d r_4_5 = R_w_4 * kinova.joint5.xyz;
    Eigen::Vector3d r_5_6 = R_w_5 * kinova.joint6.xyz;
    Eigen::Vector3d r_6_e = R_w_6 * kinova.end.xyz;

    Eigen::Vector3d r_5_e = r_5_6 + r_6_e;
    Eigen::Vector3d r_4_e = r_4_5 + r_5_e;
    Eigen::Vector3d r_3_e = r_3_4 + r_4_e;
    Eigen::Vector3d r_2_e = r_2_3 + r_3_e;
    Eigen::Vector3d r_1_e = r_1_2 + r_2_e;
    Eigen::Vector3d r_b_e = r_b_1 + r_1_e;
    Eigen::Vector3d r_w_e = r_w_b + r_b_e;

    Eigen::Vector3d P_w_1 = R_w_1 * kinova.joint1.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * kinova.joint2.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * kinova.joint3.axis;
    Eigen::Vector3d P_w_4 = R_w_4 * kinova.joint4.axis;
    Eigen::Vector3d P_w_5 = R_w_5 * kinova.joint5.axis;
    Eigen::Vector3d P_w_6 = R_w_6 * kinova.joint6.axis;

    Eigen::MatrixXd positional_jacobian(3, 6);
    positional_jacobian << skewSymmetric(P_w_1) * r_1_e,
                           skewSymmetric(P_w_2) * r_2_e,
                           skewSymmetric(P_w_3) * r_3_e,
                           skewSymmetric(P_w_4) * r_4_e,
                           skewSymmetric(P_w_5) * r_5_e,
                           skewSymmetric(P_w_6) * r_6_e;

    Eigen::MatrixXd angular_jacobian(3, 6);
    angular_jacobian << P_w_1,
                        P_w_2,
                        P_w_3,
                        P_w_4,
                        P_w_5,
                        P_w_6;

    Eigen::MatrixXd jacobian(6, 6);
    jacobian << positional_jacobian,
                angular_jacobian;

    Eigen::MatrixXd jacobian_pinv = jacobian.completeOrthogonalDecomposition().pseudoInverse();

    Eigen::Vector3d r_w_des = pos;
    Eigen::Vector3d position_error = r_w_des - r_w_e;

    Eigen::Matrix3d R_e_w = R_w_e.transpose();
    Eigen::Matrix3d R_w_des = toRot(quat[0], quat[1], quat[2], quat[3]);
    Eigen::Vector3d A_e_des = toAngleAxis(R_e_w * R_w_des);
    Eigen::Vector3d rotation_error = R_w_e * A_e_des;

    Eigen::VectorXd command_velocity_end_effector(6);
    command_velocity_end_effector << position_error,
                                     rotation_error;

    Eigen::VectorXd gv = jacobian_pinv * command_velocity_end_effector;

    return gv;
}
