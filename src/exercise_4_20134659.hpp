#pragma once

struct Kinova {
    struct Joint {
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
        Eigen::Vector3d axis;
    };

    struct Link {
        double mass;
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
        Eigen::Vector3d inertia;
    };

    Joint joint_base;
    Link link_base;
    Joint joint_1;
    Link link_1;
    Joint joint_2;
    Link link_2;
    Joint joint_3;
    Link link_3;
    Joint joint_4;
    Link link_4;
    Joint joint_5;
    Link link_5;
    Joint joint_6;
    Link link_6;
    Joint joint_end;
    Link link_end;

    Kinova() {
        joint_base.rpy = {0, 0, 0};
        joint_base.xyz = {0, 0, 0};
        joint_base.axis = {0, 0, 0};

        link_base.mass = 0.46784;
        link_base.rpy = {0, 0, 0};
        link_base.xyz = {0, 0, 0.1255};
        link_base.inertia = {0.000951270861568, 0.000951270861568, 0.000374272};

        joint_1.rpy = {0, 3.14159265359, 0};
        joint_1.xyz = {0, 0, 0.15675};
        joint_1.axis = {0, 0, 1};

        link_1.mass = 0.7477;
        link_1.rpy = {0, 0, 0};
        link_1.xyz = {0, -0.002, -0.0605};
        link_1.inertia = {0.00152031725204, 0.00152031725204, 0.00059816};

        joint_2.rpy = {-1.57079632679, 0, 3.14159265359};
        joint_2.xyz = {0, 0.0016, -0.11875};
        joint_2.axis = {0, 0, 1};

        link_2.mass = 0.99;
        link_2.rpy = {0, 0, 0};
        link_2.xyz = {0, -0.2065, -0.01};
        link_2.inertia = {0.010502207991, 0.000792, 0.010502207991};

        joint_3.rpy = {0, 3.14159265359, 0};
        joint_3.xyz = {0, -0.410, 0};
        joint_3.axis = {0, 0, 1};

        link_3.mass = 0.6763;
        link_3.rpy = {0, 0, 0};
        link_3.xyz = {0, 0.081, -0.0086};
        link_3.inertia = {0.00142022431908, 0.000304335, 0.00142022431908};

        joint_4.rpy = {-1.57079632679, 0, 3.14159265359};
        joint_4.xyz = {0, 0.2073, -0.0114};
        joint_4.axis = {0, 0, 1};

        link_4.mass = 0.463;
        link_4.rpy = {0, 0, 0};
        link_4.xyz = {0, 0.0028848942, -0.0541932613};
        link_4.inertia = {0.0004321316048, 0.0004321316048, 9.26e-05};

        joint_5.rpy = {1.57079632679, 0, 3.14159265359};
        joint_5.xyz = {0, 0, -0.10375};
        joint_5.axis = {0, 0, 1};

        link_5.mass = 0.463;
        link_5.rpy = {0, 0, 0};
        link_5.xyz = {0, 0.0497208855, -0.0028562765};
        link_5.inertia = {0.0004321316048, 9.26e-05, 0.0004321316048};

        joint_6.rpy = {-1.57079632679, 0, 3.14159265359};
        joint_6.xyz = {0, 0.10375, 0};
        joint_6.axis = {0, 0, 1};

        link_6.mass = 1.327;
        link_6.rpy = {0, 0, 0};
        link_6.xyz = {0, 0, -0.06};
        link_6.inertia = {0.0004403232387, 0.0004403232387, 0.0007416};

        joint_end.rpy = {3.14159265359, 0, 0};
        joint_end.xyz = {0, 0, -0.1600};
        joint_end.axis = {0, 0, 0};

        link_end.mass = 0.01;
        link_end.rpy = {0, 0, 0};
        link_end.xyz = {0, 0, 0};
        link_end.inertia = {0.01, 0.01, 0.01};
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

Eigen::Matrix3d toInertiaMatrix(const Eigen::Vector3d& i) {
    Eigen::Matrix3d m;
    m << i[0], 0, 0,
         0, i[1], 0,
         0, 0, i[2];
    return m;
}

// global kinova
Kinova kinova;

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
    Eigen::Matrix3d R_w_b = toRot(kinova.joint_base.rpy);
    Eigen::Matrix3d R_b_1 = toRot(kinova.joint_1.rpy);
    Eigen::Matrix3d R_1_2 = toRot(kinova.joint_2.rpy);
    Eigen::Matrix3d R_2_3 = toRot(kinova.joint_3.rpy);
    Eigen::Matrix3d R_3_4 = toRot(kinova.joint_4.rpy);
    Eigen::Matrix3d R_4_5 = toRot(kinova.joint_5.rpy);
    Eigen::Matrix3d R_5_6 = toRot(kinova.joint_6.rpy);
    Eigen::Matrix3d R_6_e = toRot(kinova.joint_end.rpy);

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

    Eigen::Vector3d r_b_b = R_w_b * kinova.link_base.xyz;
    Eigen::Vector3d r_1_1 = R_w_1 * kinova.link_1.xyz;
    Eigen::Vector3d r_2_2 = R_w_2 * kinova.link_2.xyz;
    Eigen::Vector3d r_3_3 = R_w_3 * kinova.link_3.xyz;
    Eigen::Vector3d r_4_4 = R_w_4 * kinova.link_4.xyz;
    Eigen::Vector3d r_5_5 = R_w_5 * kinova.link_5.xyz;
    Eigen::Vector3d r_6_6 = R_w_6 * kinova.link_6.xyz;
    Eigen::Vector3d r_e_e = R_w_e * kinova.link_end.xyz;

    Eigen::Vector3d r_w_b = kinova.joint_base.xyz + r_b_b;
    Eigen::Vector3d r_b_1 = R_w_b * kinova.joint_1.xyz + r_1_1;
    Eigen::Vector3d r_1_2 = R_w_1 * kinova.joint_2.xyz + r_2_2;

    Eigen::Vector3d r_2_3 = R_w_2 * kinova.joint_3.xyz + r_3_3;
    Eigen::Vector3d r_1_3 = R_w_1 * kinova.joint_2.xyz + r_2_3;

    Eigen::Vector3d r_3_4 = R_w_3 * kinova.joint_4.xyz + r_4_4;
    Eigen::Vector3d r_2_4 = R_w_2 * kinova.joint_3.xyz + r_3_4;
    Eigen::Vector3d r_1_4 = R_w_1 * kinova.joint_2.xyz + r_2_4;

    Eigen::Vector3d r_4_5 = R_w_4 * kinova.joint_5.xyz + r_5_5;
    Eigen::Vector3d r_3_5 = R_w_3 * kinova.joint_4.xyz + r_4_5;
    Eigen::Vector3d r_2_5 = R_w_2 * kinova.joint_3.xyz + r_3_5;
    Eigen::Vector3d r_1_5 = R_w_1 * kinova.joint_2.xyz + r_2_5;

    Eigen::Vector3d r_5_6 = R_w_5 * kinova.joint_6.xyz + r_6_6;
    Eigen::Vector3d r_4_6 = R_w_4 * kinova.joint_5.xyz + r_5_6;
    Eigen::Vector3d r_3_6 = R_w_3 * kinova.joint_4.xyz + r_4_6;
    Eigen::Vector3d r_2_6 = R_w_2 * kinova.joint_3.xyz + r_3_6;
    Eigen::Vector3d r_1_6 = R_w_1 * kinova.joint_2.xyz + r_2_6;

    Eigen::Vector3d r_6_e = R_w_6 * kinova.joint_end.xyz;
    Eigen::Vector3d r_5_e = R_w_5 * kinova.joint_6.xyz + r_6_e;
    Eigen::Vector3d r_4_e = R_w_4 * kinova.joint_5.xyz + r_5_e;
    Eigen::Vector3d r_3_e = R_w_3 * kinova.joint_4.xyz + r_4_e;
    Eigen::Vector3d r_2_e = R_w_2 * kinova.joint_3.xyz + r_3_e;
    Eigen::Vector3d r_1_e = R_w_1 * kinova.joint_2.xyz + r_2_e;
    Eigen::Vector3d r_b_e = R_w_b * kinova.joint_1.xyz + r_1_e;
    Eigen::Vector3d r_w_e = r_w_b + r_b_e;

    Eigen::Vector3d P_w_b = R_w_b * kinova.joint_base.axis;
    Eigen::Vector3d P_w_1 = R_w_1 * kinova.joint_1.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * kinova.joint_2.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * kinova.joint_3.axis;
    Eigen::Vector3d P_w_4 = R_w_4 * kinova.joint_4.axis;
    Eigen::Vector3d P_w_5 = R_w_5 * kinova.joint_5.axis;
    Eigen::Vector3d P_w_6 = R_w_6 * kinova.joint_6.axis;
    Eigen::Vector3d P_w_e = R_w_e * kinova.joint_end.axis;

    // positinoal jacobian
    Eigen::MatrixXd positional_jacobian_0(3, 6);
    positional_jacobian_0 << Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_1(3, 6);
    positional_jacobian_1 << skewSymmetric(P_w_1) * r_1_1,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_2(3, 6);
    positional_jacobian_2 << skewSymmetric(P_w_1) * r_1_2,
                             skewSymmetric(P_w_2) * r_2_2,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_3(3, 6);
    positional_jacobian_3 << skewSymmetric(P_w_1) * r_1_3,
                             skewSymmetric(P_w_2) * r_2_3,
                             skewSymmetric(P_w_3) * r_3_3,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_4(3, 6);
    positional_jacobian_4 << skewSymmetric(P_w_1) * r_1_4,
                             skewSymmetric(P_w_2) * r_2_4,
                             skewSymmetric(P_w_3) * r_3_4,
                             skewSymmetric(P_w_4) * r_4_4,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_5(3, 6);
    positional_jacobian_5 << skewSymmetric(P_w_1) * r_1_5,
                             skewSymmetric(P_w_2) * r_2_5,
                             skewSymmetric(P_w_3) * r_3_5,
                             skewSymmetric(P_w_4) * r_4_5,
                             skewSymmetric(P_w_5) * r_5_5,
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_6(3, 6);
    positional_jacobian_6 << skewSymmetric(P_w_1) * r_1_6,
                             skewSymmetric(P_w_2) * r_2_6,
                             skewSymmetric(P_w_3) * r_3_6,
                             skewSymmetric(P_w_4) * r_4_6,
                             skewSymmetric(P_w_5) * r_5_6,
                             skewSymmetric(P_w_6) * r_6_6;

    Eigen::MatrixXd positional_jacobian_7(3, 6);
    positional_jacobian_7 << skewSymmetric(P_w_1) * r_1_e,
                             skewSymmetric(P_w_2) * r_2_e,
                             skewSymmetric(P_w_3) * r_3_e,
                             skewSymmetric(P_w_4) * r_4_e,
                             skewSymmetric(P_w_5) * r_5_e,
                             skewSymmetric(P_w_6) * r_6_e;

    // angular jacobian
    Eigen::MatrixXd angular_jacobian_0(3, 6);
    angular_jacobian_0 << Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_1(3, 6);
    angular_jacobian_1 << P_w_1,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_2(3, 6);
    angular_jacobian_2 << P_w_1,
                          P_w_2,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_3(3, 6);
    angular_jacobian_3 << P_w_1,
                          P_w_2,
                          P_w_3,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_4(3, 6);
    angular_jacobian_4 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_5(3, 6);
    angular_jacobian_5 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          P_w_5,
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_6(3, 6);
    angular_jacobian_6 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          P_w_5,
                          P_w_6;

    Eigen::MatrixXd angular_jacobian_7(3, 6);
    angular_jacobian_7 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          P_w_5,
                          P_w_6;

    // mass matrix
    Eigen::MatrixXd mass_matrix(6, 6);
    mass_matrix.setZero();
    mass_matrix += positional_jacobian_0.transpose() * kinova.link_base.mass * positional_jacobian_0 + angular_jacobian_0.transpose() * (R_w_b * toInertiaMatrix(kinova.link_base.inertia) * R_w_b.transpose()) * angular_jacobian_0;
    mass_matrix += positional_jacobian_1.transpose() * kinova.link_1.mass * positional_jacobian_1 + angular_jacobian_1.transpose() * (R_w_1 * toInertiaMatrix(kinova.link_1.inertia) * R_w_1.transpose()) * angular_jacobian_1;
    mass_matrix += positional_jacobian_2.transpose() * kinova.link_2.mass * positional_jacobian_2 + angular_jacobian_2.transpose() * (R_w_2 * toInertiaMatrix(kinova.link_2.inertia) * R_w_2.transpose()) * angular_jacobian_2;
    mass_matrix += positional_jacobian_3.transpose() * kinova.link_3.mass * positional_jacobian_3 + angular_jacobian_3.transpose() * (R_w_3 * toInertiaMatrix(kinova.link_3.inertia) * R_w_3.transpose()) * angular_jacobian_3;
    mass_matrix += positional_jacobian_4.transpose() * kinova.link_4.mass * positional_jacobian_4 + angular_jacobian_4.transpose() * (R_w_4 * toInertiaMatrix(kinova.link_4.inertia) * R_w_4.transpose()) * angular_jacobian_4;
    mass_matrix += positional_jacobian_5.transpose() * kinova.link_5.mass * positional_jacobian_5 + angular_jacobian_5.transpose() * (R_w_5 * toInertiaMatrix(kinova.link_5.inertia) * R_w_5.transpose()) * angular_jacobian_5;
    mass_matrix += positional_jacobian_6.transpose() * kinova.link_6.mass * positional_jacobian_6 + angular_jacobian_6.transpose() * (R_w_6 * toInertiaMatrix(kinova.link_6.inertia) * R_w_6.transpose()) * angular_jacobian_6;
    mass_matrix += positional_jacobian_7.transpose() * kinova.link_end.mass * positional_jacobian_7 + angular_jacobian_7.transpose() * (R_w_e * toInertiaMatrix(kinova.link_end.inertia) * R_w_e.transpose()) * angular_jacobian_7;

    return mass_matrix;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Matrix3d R_w_b = toRot(kinova.joint_base.rpy);
    Eigen::Matrix3d R_b_1 = toRot(kinova.joint_1.rpy);
    Eigen::Matrix3d R_1_2 = toRot(kinova.joint_2.rpy);
    Eigen::Matrix3d R_2_3 = toRot(kinova.joint_3.rpy);
    Eigen::Matrix3d R_3_4 = toRot(kinova.joint_4.rpy);
    Eigen::Matrix3d R_4_5 = toRot(kinova.joint_5.rpy);
    Eigen::Matrix3d R_5_6 = toRot(kinova.joint_6.rpy);
    Eigen::Matrix3d R_6_e = toRot(kinova.joint_end.rpy);

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

    Eigen::Vector3d r_b_b = R_w_b * kinova.link_base.xyz;
    Eigen::Vector3d r_1_1 = R_w_1 * kinova.link_1.xyz;
    Eigen::Vector3d r_2_2 = R_w_2 * kinova.link_2.xyz;
    Eigen::Vector3d r_3_3 = R_w_3 * kinova.link_3.xyz;
    Eigen::Vector3d r_4_4 = R_w_4 * kinova.link_4.xyz;
    Eigen::Vector3d r_5_5 = R_w_5 * kinova.link_5.xyz;
    Eigen::Vector3d r_6_6 = R_w_6 * kinova.link_6.xyz;
    Eigen::Vector3d r_e_e = R_w_e * kinova.link_end.xyz;

    Eigen::Vector3d r_w_b = kinova.joint_base.xyz + r_b_b;
    Eigen::Vector3d r_b_1 = R_w_b * kinova.joint_1.xyz + r_1_1;
    Eigen::Vector3d r_1_2 = R_w_1 * kinova.joint_2.xyz + r_2_2;

    Eigen::Vector3d r_2_3 = R_w_2 * kinova.joint_3.xyz + r_3_3;
    Eigen::Vector3d r_1_3 = R_w_1 * kinova.joint_2.xyz + r_2_3;

    Eigen::Vector3d r_3_4 = R_w_3 * kinova.joint_4.xyz + r_4_4;
    Eigen::Vector3d r_2_4 = R_w_2 * kinova.joint_3.xyz + r_3_4;
    Eigen::Vector3d r_1_4 = R_w_1 * kinova.joint_2.xyz + r_2_4;

    Eigen::Vector3d r_4_5 = R_w_4 * kinova.joint_5.xyz + r_5_5;
    Eigen::Vector3d r_3_5 = R_w_3 * kinova.joint_4.xyz + r_4_5;
    Eigen::Vector3d r_2_5 = R_w_2 * kinova.joint_3.xyz + r_3_5;
    Eigen::Vector3d r_1_5 = R_w_1 * kinova.joint_2.xyz + r_2_5;

    Eigen::Vector3d r_5_6 = R_w_5 * kinova.joint_6.xyz + r_6_6;
    Eigen::Vector3d r_4_6 = R_w_4 * kinova.joint_5.xyz + r_5_6;
    Eigen::Vector3d r_3_6 = R_w_3 * kinova.joint_4.xyz + r_4_6;
    Eigen::Vector3d r_2_6 = R_w_2 * kinova.joint_3.xyz + r_3_6;
    Eigen::Vector3d r_1_6 = R_w_1 * kinova.joint_2.xyz + r_2_6;

    Eigen::Vector3d r_6_e = R_w_6 * kinova.joint_end.xyz;
    Eigen::Vector3d r_5_e = R_w_5 * kinova.joint_6.xyz + r_6_e;
    Eigen::Vector3d r_4_e = R_w_4 * kinova.joint_5.xyz + r_5_e;
    Eigen::Vector3d r_3_e = R_w_3 * kinova.joint_4.xyz + r_4_e;
    Eigen::Vector3d r_2_e = R_w_2 * kinova.joint_3.xyz + r_3_e;
    Eigen::Vector3d r_1_e = R_w_1 * kinova.joint_2.xyz + r_2_e;
    Eigen::Vector3d r_b_e = R_w_b * kinova.joint_1.xyz + r_1_e;
    Eigen::Vector3d r_w_e = r_w_b + r_b_e;

    Eigen::Vector3d P_w_b = R_w_b * kinova.joint_base.axis;
    Eigen::Vector3d P_w_1 = R_w_1 * kinova.joint_1.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * kinova.joint_2.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * kinova.joint_3.axis;
    Eigen::Vector3d P_w_4 = R_w_4 * kinova.joint_4.axis;
    Eigen::Vector3d P_w_5 = R_w_5 * kinova.joint_5.axis;
    Eigen::Vector3d P_w_6 = R_w_6 * kinova.joint_6.axis;
    Eigen::Vector3d P_w_e = R_w_e * kinova.joint_end.axis;

    // positinoal jacobian
    Eigen::MatrixXd positional_jacobian_0(3, 6);
    positional_jacobian_0 << Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_1(3, 6);
    positional_jacobian_1 << skewSymmetric(P_w_1) * r_1_1,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_2(3, 6);
    positional_jacobian_2 << skewSymmetric(P_w_1) * r_1_2,
                             skewSymmetric(P_w_2) * r_2_2,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_3(3, 6);
    positional_jacobian_3 << skewSymmetric(P_w_1) * r_1_3,
                             skewSymmetric(P_w_2) * r_2_3,
                             skewSymmetric(P_w_3) * r_3_3,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_4(3, 6);
    positional_jacobian_4 << skewSymmetric(P_w_1) * r_1_4,
                             skewSymmetric(P_w_2) * r_2_4,
                             skewSymmetric(P_w_3) * r_3_4,
                             skewSymmetric(P_w_4) * r_4_4,
                             Eigen::Vector3d::Zero(),
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_5(3, 6);
    positional_jacobian_5 << skewSymmetric(P_w_1) * r_1_5,
                             skewSymmetric(P_w_2) * r_2_5,
                             skewSymmetric(P_w_3) * r_3_5,
                             skewSymmetric(P_w_4) * r_4_5,
                             skewSymmetric(P_w_5) * r_5_5,
                             Eigen::Vector3d::Zero();

    Eigen::MatrixXd positional_jacobian_6(3, 6);
    positional_jacobian_6 << skewSymmetric(P_w_1) * r_1_6,
                             skewSymmetric(P_w_2) * r_2_6,
                             skewSymmetric(P_w_3) * r_3_6,
                             skewSymmetric(P_w_4) * r_4_6,
                             skewSymmetric(P_w_5) * r_5_6,
                             skewSymmetric(P_w_6) * r_6_6;

    Eigen::MatrixXd positional_jacobian_7(3, 6);
    positional_jacobian_7 << skewSymmetric(P_w_1) * r_1_e,
                             skewSymmetric(P_w_2) * r_2_e,
                             skewSymmetric(P_w_3) * r_3_e,
                             skewSymmetric(P_w_4) * r_4_e,
                             skewSymmetric(P_w_5) * r_5_e,
                             skewSymmetric(P_w_6) * r_6_e;

    // angular jacobian
    Eigen::MatrixXd angular_jacobian_0(3, 6);
    angular_jacobian_0 << Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_1(3, 6);
    angular_jacobian_1 << P_w_1,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_2(3, 6);
    angular_jacobian_2 << P_w_1,
                          P_w_2,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_3(3, 6);
    angular_jacobian_3 << P_w_1,
                          P_w_2,
                          P_w_3,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_4(3, 6);
    angular_jacobian_4 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          Eigen::Vector3d::Zero(),
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_5(3, 6);
    angular_jacobian_5 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          P_w_5,
                          Eigen::Vector3d::Zero();

    Eigen::MatrixXd angular_jacobian_6(3, 6);
    angular_jacobian_6 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          P_w_5,
                          P_w_6;

    Eigen::MatrixXd angular_jacobian_7(3, 6);
    angular_jacobian_7 << P_w_1,
                          P_w_2,
                          P_w_3,
                          P_w_4,
                          P_w_5,
                          P_w_6;

    // differential of positional_jacobian
    Eigen::Vector3d d_r_b_b = skewSymmetric(angular_jacobian_0 * gv) * R_w_b * kinova.link_base.xyz;
    Eigen::Vector3d d_r_1_1 = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.link_1.xyz;
    Eigen::Vector3d d_r_2_2 = skewSymmetric(angular_jacobian_2 * gv) * R_w_2 * kinova.link_2.xyz;
    Eigen::Vector3d d_r_3_3 = skewSymmetric(angular_jacobian_3 * gv) * R_w_3 * kinova.link_3.xyz;
    Eigen::Vector3d d_r_4_4 = skewSymmetric(angular_jacobian_4 * gv) * R_w_4 * kinova.link_4.xyz;
    Eigen::Vector3d d_r_5_5 = skewSymmetric(angular_jacobian_5 * gv) * R_w_5 * kinova.link_5.xyz;
    Eigen::Vector3d d_r_6_6 = skewSymmetric(angular_jacobian_6 * gv) * R_w_6 * kinova.link_6.xyz;
    Eigen::Vector3d d_r_e_e = skewSymmetric(angular_jacobian_7 * gv) * R_w_e * kinova.link_end.xyz;

    Eigen::Vector3d d_r_b_1 = skewSymmetric(angular_jacobian_0 * gv) * R_w_b * kinova.joint_1.xyz + d_r_1_1;
    Eigen::Vector3d d_r_1_2 = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.joint_2.xyz + d_r_2_2;

    Eigen::Vector3d d_r_2_3 = skewSymmetric(angular_jacobian_2 * gv) * R_w_2 * kinova.joint_3.xyz + d_r_3_3;
    Eigen::Vector3d d_r_1_3 = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.joint_2.xyz + d_r_2_3;

    Eigen::Vector3d d_r_3_4 = skewSymmetric(angular_jacobian_3 * gv) * R_w_3 * kinova.joint_4.xyz + d_r_4_4;
    Eigen::Vector3d d_r_2_4 = skewSymmetric(angular_jacobian_2 * gv) * R_w_2 * kinova.joint_3.xyz + d_r_3_4;
    Eigen::Vector3d d_r_1_4 = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.joint_2.xyz + d_r_2_4;

    Eigen::Vector3d d_r_4_5 = skewSymmetric(angular_jacobian_4 * gv) * R_w_4 * kinova.joint_5.xyz + d_r_5_5;
    Eigen::Vector3d d_r_3_5 = skewSymmetric(angular_jacobian_3 * gv) * R_w_3 * kinova.joint_4.xyz + d_r_4_5;
    Eigen::Vector3d d_r_2_5 = skewSymmetric(angular_jacobian_2 * gv) * R_w_2 * kinova.joint_3.xyz + d_r_3_5;
    Eigen::Vector3d d_r_1_5 = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.joint_2.xyz + d_r_2_5;

    Eigen::Vector3d d_r_5_6 = skewSymmetric(angular_jacobian_5 * gv) * R_w_5 * kinova.joint_6.xyz + d_r_6_6;
    Eigen::Vector3d d_r_4_6 = skewSymmetric(angular_jacobian_4 * gv) * R_w_4 * kinova.joint_5.xyz + d_r_5_6;
    Eigen::Vector3d d_r_3_6 = skewSymmetric(angular_jacobian_3 * gv) * R_w_3 * kinova.joint_4.xyz + d_r_4_6;
    Eigen::Vector3d d_r_2_6 = skewSymmetric(angular_jacobian_2 * gv) * R_w_2 * kinova.joint_3.xyz + d_r_3_6;
    Eigen::Vector3d d_r_1_6 = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.joint_2.xyz + d_r_2_6;

    Eigen::Vector3d d_r_6_e = skewSymmetric(angular_jacobian_6 * gv) * R_w_6 * kinova.joint_end.xyz;
    Eigen::Vector3d d_r_5_e = skewSymmetric(angular_jacobian_5 * gv) * R_w_5 * kinova.joint_6.xyz + d_r_6_e;
    Eigen::Vector3d d_r_4_e = skewSymmetric(angular_jacobian_4 * gv) * R_w_4 * kinova.joint_5.xyz + d_r_5_e;
    Eigen::Vector3d d_r_3_e = skewSymmetric(angular_jacobian_3 * gv) * R_w_3 * kinova.joint_4.xyz + d_r_4_e;
    Eigen::Vector3d d_r_2_e = skewSymmetric(angular_jacobian_2 * gv) * R_w_2 * kinova.joint_3.xyz + d_r_3_e;
    Eigen::Vector3d d_r_1_e = skewSymmetric(angular_jacobian_1 * gv) * R_w_1 * kinova.joint_2.xyz + d_r_2_e;
    Eigen::Vector3d d_r_b_e = skewSymmetric(angular_jacobian_0 * gv) * R_w_b * kinova.joint_1.xyz + d_r_1_e;

    Eigen::MatrixXd d_positional_jacobian_0(3, 6);
    d_positional_jacobian_0 << Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_positional_jacobian_1(3, 6);
    d_positional_jacobian_1 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_1 + skewSymmetric(P_w_1) * d_r_1_1,
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_positional_jacobian_2(3, 6);
    d_positional_jacobian_2 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_2 + skewSymmetric(P_w_1) * d_r_1_2,
                               skewSymmetric(skewSymmetric(angular_jacobian_2 * gv) * P_w_2) * r_2_2 + skewSymmetric(P_w_2) * d_r_2_2,
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_positional_jacobian_3(3, 6);
    d_positional_jacobian_3 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_3 + skewSymmetric(P_w_1) * d_r_1_3,
                               skewSymmetric(skewSymmetric(angular_jacobian_2 * gv) * P_w_2) * r_2_3 + skewSymmetric(P_w_2) * d_r_2_3,
                               skewSymmetric(skewSymmetric(angular_jacobian_3 * gv) * P_w_3) * r_3_3 + skewSymmetric(P_w_3) * d_r_3_3,
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_positional_jacobian_4(3, 6);
    d_positional_jacobian_4 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_4 + skewSymmetric(P_w_1) * d_r_1_4,
                               skewSymmetric(skewSymmetric(angular_jacobian_2 * gv) * P_w_2) * r_2_4 + skewSymmetric(P_w_2) * d_r_2_4,
                               skewSymmetric(skewSymmetric(angular_jacobian_3 * gv) * P_w_3) * r_3_4 + skewSymmetric(P_w_3) * d_r_3_4,
                               skewSymmetric(skewSymmetric(angular_jacobian_4 * gv) * P_w_4) * r_4_4 + skewSymmetric(P_w_4) * d_r_4_4,
                               Eigen::Vector3d::Zero(),
                               Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_positional_jacobian_5(3, 6);
    d_positional_jacobian_5 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_5 + skewSymmetric(P_w_1) * d_r_1_5,
                               skewSymmetric(skewSymmetric(angular_jacobian_2 * gv) * P_w_2) * r_2_5 + skewSymmetric(P_w_2) * d_r_2_5,
                               skewSymmetric(skewSymmetric(angular_jacobian_3 * gv) * P_w_3) * r_3_5 + skewSymmetric(P_w_3) * d_r_3_5,
                               skewSymmetric(skewSymmetric(angular_jacobian_4 * gv) * P_w_4) * r_4_5 + skewSymmetric(P_w_4) * d_r_4_5,
                               skewSymmetric(skewSymmetric(angular_jacobian_5 * gv) * P_w_5) * r_5_5 + skewSymmetric(P_w_5) * d_r_5_5,
                               Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_positional_jacobian_6(3, 6);
    d_positional_jacobian_6 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_6 + skewSymmetric(P_w_1) * d_r_1_6,
                               skewSymmetric(skewSymmetric(angular_jacobian_2 * gv) * P_w_2) * r_2_6 + skewSymmetric(P_w_2) * d_r_2_6,
                               skewSymmetric(skewSymmetric(angular_jacobian_3 * gv) * P_w_3) * r_3_6 + skewSymmetric(P_w_3) * d_r_3_6,
                               skewSymmetric(skewSymmetric(angular_jacobian_4 * gv) * P_w_4) * r_4_6 + skewSymmetric(P_w_4) * d_r_4_6,
                               skewSymmetric(skewSymmetric(angular_jacobian_5 * gv) * P_w_5) * r_5_6 + skewSymmetric(P_w_5) * d_r_5_6,
                               skewSymmetric(skewSymmetric(angular_jacobian_6 * gv) * P_w_6) * r_6_6 + skewSymmetric(P_w_6) * d_r_6_6;

    Eigen::MatrixXd d_positional_jacobian_7(3, 6);
    d_positional_jacobian_7 << skewSymmetric(skewSymmetric(angular_jacobian_1 * gv) * P_w_1) * r_1_e + skewSymmetric(P_w_1) * d_r_1_e,
                               skewSymmetric(skewSymmetric(angular_jacobian_2 * gv) * P_w_2) * r_2_e + skewSymmetric(P_w_2) * d_r_2_e,
                               skewSymmetric(skewSymmetric(angular_jacobian_3 * gv) * P_w_3) * r_3_e + skewSymmetric(P_w_3) * d_r_3_e,
                               skewSymmetric(skewSymmetric(angular_jacobian_4 * gv) * P_w_4) * r_4_e + skewSymmetric(P_w_4) * d_r_4_e,
                               skewSymmetric(skewSymmetric(angular_jacobian_5 * gv) * P_w_5) * r_5_e + skewSymmetric(P_w_5) * d_r_5_e,
                               skewSymmetric(skewSymmetric(angular_jacobian_6 * gv) * P_w_6) * r_6_e + skewSymmetric(P_w_6) * d_r_6_e;

    // differential of angular jacobian
    Eigen::MatrixXd d_angular_jacobian_0(3, 6);
    d_angular_jacobian_0 << Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_angular_jacobian_1(3, 6);
    d_angular_jacobian_1 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_angular_jacobian_2(3, 6);
    d_angular_jacobian_2 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            skewSymmetric(angular_jacobian_2 * gv) * P_w_2,
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_angular_jacobian_3(3, 6);
    d_angular_jacobian_3 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            skewSymmetric(angular_jacobian_2 * gv) * P_w_2,
                            skewSymmetric(angular_jacobian_3 * gv) * P_w_3,
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_angular_jacobian_4(3, 6);
    d_angular_jacobian_4 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            skewSymmetric(angular_jacobian_2 * gv) * P_w_2,
                            skewSymmetric(angular_jacobian_3 * gv) * P_w_3,
                            skewSymmetric(angular_jacobian_4 * gv) * P_w_4,
                            Eigen::Vector3d::Zero(),
                            Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_angular_jacobian_5(3, 6);
    d_angular_jacobian_5 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            skewSymmetric(angular_jacobian_2 * gv) * P_w_2,
                            skewSymmetric(angular_jacobian_3 * gv) * P_w_3,
                            skewSymmetric(angular_jacobian_4 * gv) * P_w_4,
                            skewSymmetric(angular_jacobian_5 * gv) * P_w_5,
                            Eigen::Vector3d::Zero();

    Eigen::MatrixXd d_angular_jacobian_6(3, 6);
    d_angular_jacobian_6 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            skewSymmetric(angular_jacobian_2 * gv) * P_w_2,
                            skewSymmetric(angular_jacobian_3 * gv) * P_w_3,
                            skewSymmetric(angular_jacobian_4 * gv) * P_w_4,
                            skewSymmetric(angular_jacobian_5 * gv) * P_w_5,
                            skewSymmetric(angular_jacobian_6 * gv) * P_w_6;

    Eigen::MatrixXd d_angular_jacobian_7(3, 6);
    d_angular_jacobian_7 << skewSymmetric(angular_jacobian_1 * gv) * P_w_1,
                            skewSymmetric(angular_jacobian_2 * gv) * P_w_2,
                            skewSymmetric(angular_jacobian_3 * gv) * P_w_3,
                            skewSymmetric(angular_jacobian_4 * gv) * P_w_4,
                            skewSymmetric(angular_jacobian_5 * gv) * P_w_5,
                            skewSymmetric(angular_jacobian_6 * gv) * P_w_6;

    // nonlinearities
    Eigen::Vector3d gravity = {0, 0, -9.81};
    Eigen::VectorXd nonlinearities(6);
    nonlinearities.setZero();
    nonlinearities += positional_jacobian_0.transpose() * kinova.link_base.mass * d_positional_jacobian_0 * gv
                    + angular_jacobian_0.transpose() * (R_w_b * toInertiaMatrix(kinova.link_base.inertia) * R_w_b.transpose()) * d_angular_jacobian_0 * gv
                    + angular_jacobian_0.transpose() * skewSymmetric(angular_jacobian_0 * gv) * (R_w_b * toInertiaMatrix(kinova.link_base.inertia) * R_w_b.transpose()) * (angular_jacobian_0 * gv)
                    - positional_jacobian_0.transpose() * kinova.link_base.mass * gravity;
    nonlinearities += positional_jacobian_1.transpose() * kinova.link_1.mass * d_positional_jacobian_1 * gv
                    + angular_jacobian_1.transpose() * (R_w_1 * toInertiaMatrix(kinova.link_1.inertia) * R_w_1.transpose()) * d_angular_jacobian_1 * gv
                    + angular_jacobian_1.transpose() * skewSymmetric(angular_jacobian_1 * gv) * (R_w_1 * toInertiaMatrix(kinova.link_1.inertia) * R_w_1.transpose()) * (angular_jacobian_1 * gv)
                    - positional_jacobian_1.transpose() * kinova.link_1.mass * gravity;
    nonlinearities += positional_jacobian_2.transpose() * kinova.link_2.mass * d_positional_jacobian_2 * gv
                    + angular_jacobian_2.transpose() * (R_w_2 * toInertiaMatrix(kinova.link_2.inertia) * R_w_2.transpose()) * d_angular_jacobian_2 * gv
                    + angular_jacobian_2.transpose() * skewSymmetric(angular_jacobian_2 * gv) * (R_w_2 * toInertiaMatrix(kinova.link_2.inertia) * R_w_2.transpose()) * (angular_jacobian_2 * gv)
                    - positional_jacobian_2.transpose() * kinova.link_2.mass * gravity;
    nonlinearities += positional_jacobian_3.transpose() * kinova.link_3.mass * d_positional_jacobian_3 * gv
                    + angular_jacobian_3.transpose() * (R_w_3 * toInertiaMatrix(kinova.link_3.inertia) * R_w_3.transpose()) * d_angular_jacobian_3 * gv
                    + angular_jacobian_3.transpose() * skewSymmetric(angular_jacobian_3 * gv) * (R_w_3 * toInertiaMatrix(kinova.link_3.inertia) * R_w_3.transpose()) * (angular_jacobian_3 * gv)
                    - positional_jacobian_3.transpose() * kinova.link_3.mass * gravity;
    nonlinearities += positional_jacobian_4.transpose() * kinova.link_4.mass * d_positional_jacobian_4 * gv
                    + angular_jacobian_4.transpose() * (R_w_4 * toInertiaMatrix(kinova.link_4.inertia) * R_w_4.transpose()) * d_angular_jacobian_4 * gv
                    + angular_jacobian_4.transpose() * skewSymmetric(angular_jacobian_4 * gv) * (R_w_4 * toInertiaMatrix(kinova.link_4.inertia) * R_w_4.transpose()) * (angular_jacobian_4 * gv)
                    - positional_jacobian_4.transpose() * kinova.link_4.mass * gravity;
    nonlinearities += positional_jacobian_5.transpose() * kinova.link_5.mass * d_positional_jacobian_5 * gv
                    + angular_jacobian_5.transpose() * (R_w_5 * toInertiaMatrix(kinova.link_5.inertia) * R_w_5.transpose()) * d_angular_jacobian_5 * gv
                    + angular_jacobian_5.transpose() * skewSymmetric(angular_jacobian_5 * gv) * (R_w_5 * toInertiaMatrix(kinova.link_5.inertia) * R_w_5.transpose()) * (angular_jacobian_5 * gv)
                    - positional_jacobian_5.transpose() * kinova.link_5.mass * gravity;
    nonlinearities += positional_jacobian_6.transpose() * kinova.link_6.mass * d_positional_jacobian_6 * gv
                    + angular_jacobian_6.transpose() * (R_w_6 * toInertiaMatrix(kinova.link_6.inertia) * R_w_6.transpose()) * d_angular_jacobian_6 * gv
                    + angular_jacobian_6.transpose() * skewSymmetric(angular_jacobian_6 * gv) * (R_w_6 * toInertiaMatrix(kinova.link_6.inertia) * R_w_6.transpose()) * (angular_jacobian_6 * gv)
                    - positional_jacobian_6.transpose() * kinova.link_6.mass * gravity;
    nonlinearities += positional_jacobian_7.transpose() * kinova.link_end.mass * d_positional_jacobian_7 * gv
                    + angular_jacobian_7.transpose() * (R_w_e * toInertiaMatrix(kinova.link_end.inertia) * R_w_e.transpose()) * d_angular_jacobian_7 * gv
                    + angular_jacobian_7.transpose() * skewSymmetric(angular_jacobian_7 * gv) * (R_w_e * toInertiaMatrix(kinova.link_end.inertia) * R_w_e.transpose()) * (angular_jacobian_7 * gv)
                    - positional_jacobian_7.transpose() * kinova.link_end.mass * gravity;

    return nonlinearities;
}

