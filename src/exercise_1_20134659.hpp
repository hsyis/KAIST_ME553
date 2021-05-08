#pragma once

struct Kinova {
    struct Joint {
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
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

        joint1.rpy = {0, 3.14159265359, 0};
        joint1.xyz = {0, 0, 0.15675};

        joint2.rpy = {-1.57079632679, 0, 3.14159265359};
        joint2.xyz = {0, 0.0016, -0.11875};

        joint3.rpy = {0, 3.14159265359, 0};
        joint3.xyz = {0, -0.410, 0};

        joint4.rpy = {-1.57079632679, 0, 3.14159265359};
        joint4.xyz = {0, 0.2073, -0.0114};

        joint5.rpy = {1.57079632679, 0, 3.14159265359};
        joint5.xyz = {0, 0, -0.10375};

        joint6.rpy = {-1.57079632679, 0, 3.14159265359};
        joint6.xyz = {0, 0.10375, 0};

        end.rpy = {3.14159265359, 0, 0};
        end.xyz = {0, 0, -0.1600};
    }
};

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

Eigen::Matrix3d toRotYaw(double theta) {
    Eigen::VectorXd rpy(3);
    rpy << 0, 0, theta;
    return toRot(rpy);
}

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
    Kinova kinova;

    Eigen::Matrix3d R_w_b = toRot(kinova.base.rpy);
    Eigen::Matrix3d R_b_1 = toRot(kinova.joint1.rpy);
    Eigen::Matrix3d R_1_2 = toRot(kinova.joint2.rpy);
    Eigen::Matrix3d R_2_3 = toRot(kinova.joint3.rpy);
    Eigen::Matrix3d R_3_4 = toRot(kinova.joint4.rpy);
    Eigen::Matrix3d R_4_5 = toRot(kinova.joint5.rpy);
    Eigen::Matrix3d R_5_6 = toRot(kinova.joint6.rpy);
    Eigen::Matrix3d R_6_e = toRot(kinova.end.rpy);

    Eigen::Matrix3d R_1_1 = toRotYaw(gc[0]);
    Eigen::Matrix3d R_2_2 = toRotYaw(gc[1]);
    Eigen::Matrix3d R_3_3 = toRotYaw(gc[2]);
    Eigen::Matrix3d R_4_4 = toRotYaw(gc[3]);
    Eigen::Matrix3d R_5_5 = toRotYaw(gc[4]);
    Eigen::Matrix3d R_6_6 = toRotYaw(gc[5]);

    Eigen::Matrix3d R_w_1 = R_w_b * R_b_1 * R_1_1;
    Eigen::Matrix3d R_w_2 = R_w_1 * R_1_2 * R_2_2;
    Eigen::Matrix3d R_w_3 = R_w_2 * R_2_3 * R_3_3;
    Eigen::Matrix3d R_w_4 = R_w_3 * R_3_4 * R_4_4;
    Eigen::Matrix3d R_w_5 = R_w_4 * R_4_5 * R_5_5;
    Eigen::Matrix3d R_w_6 = R_w_5 * R_5_6 * R_6_6;

    Eigen::Vector3d r_b_1 = R_w_b * kinova.joint1.xyz;
    Eigen::Vector3d r_1_2 = R_w_1 * kinova.joint2.xyz;
    Eigen::Vector3d r_2_3 = R_w_2 * kinova.joint3.xyz;
    Eigen::Vector3d r_3_4 = R_w_3 * kinova.joint4.xyz;
    Eigen::Vector3d r_4_5 = R_w_4 * kinova.joint5.xyz;
    Eigen::Vector3d r_5_6 = R_w_5 * kinova.joint6.xyz;
    Eigen::Vector3d r_6_e = R_w_6 * kinova.end.xyz;

    return r_b_1 + r_1_2 + r_2_3 + r_3_4 + r_4_5 + r_5_6 + r_6_e;
}
