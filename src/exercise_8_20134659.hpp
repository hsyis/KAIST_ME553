#pragma once

struct Anymal {
    struct Joint {
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
        Eigen::Vector3d axis;
    };

    struct Link {
        double mass;
        Eigen::Vector3d rpy;
        Eigen::Vector3d xyz;
        Eigen::Matrix<double, 6, 1> inertia;
    };

    Joint joint_base;
    Link link_base;
    Joint joint_hip;
    Link link_hip;
    Joint joint_thigh;
    Link link_thigh;
    Joint joint_shank;
    Link link_shank;
    Joint joint_adapter; // shank_to_adapter
    Link link_adapter;
    Joint joint_foot;    // adapter_to_foot

    Anymal() {
        joint_base.rpy = {0.0, 0.0, 0.0};
        joint_base.xyz = {0.0, 0.0, 0.0};
        joint_base.axis = {0, 0, 0};

        link_base.mass = 16.793507758;
        link_base.rpy = {0.0, 0.0, 0.0};
        link_base.xyz = {-0.001960558279, -0.001413217745, 0.050207125344};
        link_base.inertia << 0.217391101503, -0.00132873239126, -0.00228200226173, 0.639432546734, -0.00138078263145, 0.62414077654;

        joint_hip.rpy = {0.0, 0.0, 0.0};
        joint_hip.xyz = {0.277, 0.116, 0.0};
        joint_hip.axis = {1, 0, 0};

        link_hip.mass = 1.42462064;
        link_hip.rpy = {0.0, 0.0, 0.0};
        link_hip.xyz = {0.064516258147, -0.003787101702, -0.000152184388};
        link_hip.inertia << 0.00243023349564, -1.53023971e-05, -2.1819095354e-05, 0.00230257239103, 2.6473021273e-05, 0.0019806759227;

        joint_thigh.rpy = {0.0, 0.0, 0.0};
        joint_thigh.xyz = {0.0635, 0.041, 0.0};
        joint_thigh.axis = {0, 1, 0};

        link_thigh.mass = 1.634976467;
        link_thigh.rpy = {0.0, 0.0, 0.0};
        link_thigh.xyz = {-0.003897968082, 0.054226618537, -0.214583373795};
        link_thigh.inertia << 0.0120367944369, 6.762065206e-05, 0.000287806340448, 0.0120643637939, -0.00140610131218, 0.00249422574881;

        joint_shank.rpy = {0.0, 0.0, 0.0};
        joint_shank.xyz = {0.0, 0.109, -0.25};
        joint_shank.axis = {0, 1, 0};

        link_shank.mass = 0.207204302;
        link_shank.rpy = {0.0, 0.0, 0.0};
        link_shank.xyz = {0.030816858139, -0.004617229294, 0.000893125713};
        link_shank.inertia << 0.0002104880248, -5.6750980345e-05, 1.0127699391e-05, 0.000676270210023, -8.22869024e-07, 0.000545032674924;

        joint_adapter.rpy = {0.0, 0.0, 0.0};
        joint_adapter.xyz = {0.1, -0.02, 0.0};
        joint_adapter.axis = {0, 0, 0};

        link_adapter.mass = 0.140170767;
        link_adapter.rpy = {0.0, 0.0, 0.0};
        link_adapter.xyz = {-8.66e-10, -1.472e-09, -0.244345749188};
        link_adapter.inertia << 0.00159938741862, -9.32e-13, 1.039e-11, 0.00159938741932, 1.7563e-11, 5.4423177329e-05;

        joint_foot.rpy = {0.0, 0.0, 0.0};
        joint_foot.xyz = {0.0, 0.0, -0.32125};
        joint_foot.axis = {0, 0, 0};
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

Eigen::Matrix3d toInertiaMatrix(const Eigen::VectorXd& i) {
    Eigen::Matrix3d m;

    switch (i.size()) {
    case 3:
        m << i[0], 0, 0,
             0, i[1], 0,
             0, 0, i[2];
        break;
    case 6:
        m << i[0], i[1], i[2],
             i[1], i[3], i[4],
             i[2], i[4], i[5];
        break;
    }

    return m;
}

Eigen::Vector3d getCompPos(double m1, double m2, const Eigen::Vector3d& r1, const Eigen::Vector3d& r2) {
    return (m1 * r1 + m2 * r2) / (m1 + m2);
}

Eigen::Matrix3d getCompInertia(double m1, double m2, const Eigen::Matrix3d& i1, const Eigen::Matrix3d& i2, const Eigen::Vector3d& r1, const Eigen::Vector3d& r2) {
    return i1 + i2 - (m1 * skewSymmetric(r1) * skewSymmetric(r1)) - (m2 * skewSymmetric(r2) * skewSymmetric(r2));
}

// global anymal
Anymal anymal;

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrixUsingCRBA (const Eigen::VectorXd& gc) {
    /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

    Eigen::Matrix3d R_w_b = toRot(gc[3], gc[4], gc[5], gc[6]);
    Eigen::Matrix3d R_b_1 = toRot(anymal.joint_hip.rpy);
    Eigen::Matrix3d R_1_2 = toRot(anymal.joint_thigh.rpy);
    Eigen::Matrix3d R_2_3 = toRot(anymal.joint_shank.rpy);
    Eigen::Matrix3d R_3_4 = toRot(anymal.joint_adapter.rpy);

    Eigen::Matrix3d R_1_1 = toRot(gc[7], Rotation::roll);
    Eigen::Matrix3d R_2_2 = toRot(gc[8], Rotation::pitch);
    Eigen::Matrix3d R_3_3 = toRot(gc[9], Rotation::pitch);

    Eigen::Matrix3d R_w_1 = R_w_b * R_b_1 * R_1_1;
    Eigen::Matrix3d R_w_2 = R_w_1 * R_1_2 * R_2_2;
    Eigen::Matrix3d R_w_3 = R_w_2 * R_2_3 * R_3_3;
    Eigen::Matrix3d R_w_4 = R_w_3 * R_3_4;

    Eigen::Vector3d r_b_bCOM = R_w_b * anymal.link_base.xyz;
    Eigen::Vector3d r_1_1COM = R_w_1 * anymal.link_hip.xyz;
    Eigen::Vector3d r_2_2COM = R_w_2 * anymal.link_thigh.xyz;
    Eigen::Vector3d r_3_3COM = R_w_3 * anymal.link_shank.xyz;
    Eigen::Vector3d r_4_4COM = R_w_4 * anymal.link_adapter.xyz;

    Eigen::Vector3d r_w_b = {gc[0], gc[1], gc[2]};
    Eigen::Vector3d r_b_1 = R_w_b * anymal.joint_hip.xyz;
    Eigen::Vector3d r_1_2 = R_w_1 * anymal.joint_thigh.xyz;
    Eigen::Vector3d r_2_3 = R_w_2 * anymal.joint_shank.xyz;
    Eigen::Vector3d r_3_4 = R_w_3 * anymal.joint_adapter.xyz;
    Eigen::Vector3d r_4_e = R_w_4 * anymal.joint_foot.xyz;

    Eigen::Vector3d r_w_bCOM = r_w_b + r_b_bCOM;
    Eigen::Vector3d r_b_1COM = r_b_1 + r_1_1COM;
    Eigen::Vector3d r_1_2COM = r_1_2 + r_2_2COM;
    Eigen::Vector3d r_2_3COM = r_2_3 + r_3_3COM;
    Eigen::Vector3d r_3_4COM = r_3_4 + r_4_4COM;

    Eigen::Vector3d r_w_1 = r_w_b + r_b_1;
    Eigen::Vector3d r_w_2 = r_w_1 + r_1_2;
    Eigen::Vector3d r_w_3 = r_w_2 + r_2_3;
    Eigen::Vector3d r_w_4 = r_w_3 + r_3_4;
    Eigen::Vector3d r_w_e = r_w_4 + r_4_e;

    Eigen::Vector3d r_w_1COM = r_w_b + r_b_1COM;
    Eigen::Vector3d r_w_2COM = r_w_1 + r_1_2COM;
    Eigen::Vector3d r_w_3COM = r_w_2 + r_2_3COM;
    Eigen::Vector3d r_w_4COM = r_w_3 + r_3_4COM;

    Eigen::Vector3d r_3_e = r_3_4 + r_4_e;
    Eigen::Vector3d r_2_e = r_2_3 + r_3_e;
    Eigen::Vector3d r_1_e = r_1_2 + r_2_e;
    Eigen::Vector3d r_b_e = r_b_1 + r_1_e;

    Eigen::Vector3d P_w_b = R_w_b * anymal.joint_base.axis;
    Eigen::Vector3d P_w_1 = R_w_1 * anymal.joint_hip.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * anymal.joint_thigh.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * anymal.joint_shank.axis;

    Eigen::Matrix3d I_b = R_w_b * toInertiaMatrix(anymal.link_base.inertia) * R_w_b.transpose();
    Eigen::Matrix3d I_1 = R_w_1 * toInertiaMatrix(anymal.link_hip.inertia) * R_w_1.transpose();
    Eigen::Matrix3d I_2 = R_w_2 * toInertiaMatrix(anymal.link_thigh.inertia) * R_w_2.transpose();
    Eigen::Matrix3d I_3 = R_w_3 * toInertiaMatrix(anymal.link_shank.inertia) * R_w_3.transpose();
    Eigen::Matrix3d I_4 = R_w_4 * toInertiaMatrix(anymal.link_adapter.inertia) * R_w_4.transpose();

    double M_b = anymal.link_base.mass;
    double M_1 = anymal.link_hip.mass;
    double M_2 = anymal.link_thigh.mass;
    double M_3 = anymal.link_shank.mass;
    double M_4 = anymal.link_adapter.mass;

    Eigen::MatrixXd spatial_inertia_matrix(6, 6);
    {
        spatial_inertia_matrix.setZero();

        Eigen::Vector3d r_w_b1 = getCompPos(M_b, M_1, r_w_bCOM, r_w_1COM);
        Eigen::Matrix3d I_b1 = getCompInertia(M_b, M_1, I_b, I_1, r_w_bCOM - r_w_b1, r_w_1COM - r_w_b1);
        double M_b1 = M_b + M_1;

        Eigen::Vector3d r_w_b12 = getCompPos(M_b1, M_2, r_w_b1, r_w_2COM);
        Eigen::Matrix3d I_b12 = getCompInertia(M_b1, M_2, I_b1, I_2, r_w_b1 - r_w_b12, r_w_2COM - r_w_b12);
        double M_b12 = M_b + M_1 + M_2;

        Eigen::Vector3d r_w_b123 = getCompPos(M_b12, M_3, r_w_b12, r_w_3COM);
        Eigen::Matrix3d I_b123 = getCompInertia(M_b12, M_3, I_b12, I_3, r_w_b12 - r_w_b123, r_w_3COM - r_w_b123);
        double M_b123 = M_b + M_1 + M_2 + M_3;

        Eigen::Vector3d r_w_b1234 = getCompPos(M_b123, M_4, r_w_b123, r_w_4COM);
        Eigen::Matrix3d I_b1234 = getCompInertia(M_b123, M_4, I_b123, I_4, r_w_b123 - r_w_b1234, r_w_4COM - r_w_b1234);
        double M_b1234 = M_b + M_1 + M_2 + M_3 + M_4;

        spatial_inertia_matrix << M_b1234 * Eigen::Matrix3d::Identity(), -M_b1234 * skewSymmetric(r_w_b1234 - r_w_b),
                                  M_b1234 * skewSymmetric(r_w_b1234 - r_w_b), I_b1234 - M_b1234 * skewSymmetric(r_w_b1234 - r_w_b) * skewSymmetric(r_w_b1234 - r_w_b);
    }

    Eigen::MatrixXd mass_matrix_1(3, 6);
    {
        mass_matrix_1.setZero();

        //1
        Eigen::Vector3d r_w_12 = getCompPos(M_1, M_2, r_w_1COM, r_w_2COM);
        Eigen::Matrix3d I_12 = getCompInertia(M_1, M_2, I_1, I_2, r_w_1COM - r_w_12, r_w_2COM - r_w_12);
        double M_12 = M_1 + M_2;

        Eigen::Vector3d r_w_34 = getCompPos(M_3, M_4, r_w_3COM, r_w_4COM);
        Eigen::Matrix3d I_34 = getCompInertia(M_3, M_4, I_3, I_4, r_w_3COM - r_w_34, r_w_4COM - r_w_34);
        double M_34 = M_3 + M_4;

        Eigen::Vector3d r_w_1234 = getCompPos(M_12, M_34, r_w_12, r_w_34);
        Eigen::Matrix3d I_1234 = getCompInertia(M_12, M_34, I_12, I_34, r_w_12 - r_w_1234, r_w_34 - r_w_1234);
        double M_1234 = M_12 + M_34;

        Eigen::MatrixXd spatial_inertia_1(6, 6);
        spatial_inertia_1 << M_1234 * Eigen::Matrix3d::Identity(), -M_1234 * skewSymmetric(r_w_1234 - r_w_1),
                             M_1234 * skewSymmetric(r_w_1234 - r_w_1), I_1234 - M_1234 * skewSymmetric(r_w_1234 - r_w_1) * skewSymmetric(r_w_1234 - r_w_1);

        Eigen::VectorXd fictitious_force_1(6);
        fictitious_force_1 << M_1234 * skewSymmetric(Eigen::Vector3d::Zero()) * skewSymmetric(Eigen::Vector3d::Zero()) * (r_w_1234 - r_w_1),
                              skewSymmetric(Eigen::Vector3d::Zero()) * (I_1234 - M_1234 * skewSymmetric(r_w_1234 - r_w_1) * skewSymmetric(r_w_1234 - r_w_1)) * Eigen::Vector3d::Zero();

        Eigen::VectorXd u_dot_1(6);
        u_dot_1 << Eigen::Vector3d::UnitX(),
                   Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_1 = spatial_inertia_1 * u_dot_1 + fictitious_force_1;

        Eigen::VectorXd s_1(6);
        s_1 << Eigen::Vector3d::Zero(),
               P_w_1;

        double mm1 = (s_1.transpose() * tau_1)(0);

        //2
        Eigen::VectorXd u_dot_2(6);
        u_dot_2 << Eigen::Vector3d::UnitY(),
                   Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_2 = spatial_inertia_1 * u_dot_2 + fictitious_force_1;

        Eigen::VectorXd s_2(6);
        s_2 << Eigen::Vector3d::Zero(),
               P_w_1;

        double mm2 = (s_2.transpose() * tau_2)(0);

        //3
        Eigen::VectorXd u_dot_3(6);
        u_dot_3 << Eigen::Vector3d::UnitZ(),
                   Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_3 = spatial_inertia_1 * u_dot_3 + fictitious_force_1;

        Eigen::VectorXd s_3(6);
        s_3 << Eigen::Vector3d::Zero(),
               P_w_1;

        double mm3 = (s_3.transpose() * tau_3)(0);

        //4
        Eigen::VectorXd u_dot_4(6);
        u_dot_4 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitX()) * (r_w_1 - r_w_b),
                   Eigen::Vector3d::UnitX();

        Eigen::VectorXd tau_4 = spatial_inertia_1 * u_dot_4 + fictitious_force_1;

        Eigen::VectorXd s_4(6);
        s_4 << Eigen::Vector3d::Zero(),
               P_w_1;

        double mm4 = (s_4.transpose() * tau_4)(0);

        //5
        Eigen::VectorXd u_dot_5(6);
        u_dot_5 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitY()) * (r_w_1 - r_w_b),
                   Eigen::Vector3d::UnitY();

        Eigen::VectorXd tau_5 = spatial_inertia_1 * u_dot_5 + fictitious_force_1;

        Eigen::VectorXd s_5(6);
        s_5 << Eigen::Vector3d::Zero(),
               P_w_1;

        double mm5 = (s_5.transpose() * tau_5)(0);

        //6
        Eigen::VectorXd u_dot_6(6);
        u_dot_6 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitZ()) * (r_w_1 - r_w_b),
                   Eigen::Vector3d::UnitZ();

        Eigen::VectorXd tau_6 = spatial_inertia_1 * u_dot_6 + fictitious_force_1;

        Eigen::VectorXd s_6(6);
        s_6 << Eigen::Vector3d::Zero(),
               P_w_1;

        double mm6 = (s_6.transpose() * tau_6)(0);

        //7
        Eigen::Vector3d r_w_234 = getCompPos(M_2, M_34, r_w_2COM, r_w_34);
        Eigen::Matrix3d I_234 = getCompInertia(M_2, M_34, I_2, I_34, r_w_2COM - r_w_234, r_w_34 - r_w_234);
        double M_234 = M_2 + M_34;

        Eigen::MatrixXd spatial_inertia_2(6, 6);
        spatial_inertia_2 << M_234 * Eigen::Matrix3d::Identity(), -M_234 * skewSymmetric(r_w_234 - r_w_2),
                             M_234 * skewSymmetric(r_w_234 - r_w_2), I_234 - M_234 * skewSymmetric(r_w_234 - r_w_2) * skewSymmetric(r_w_234 - r_w_2);

        Eigen::VectorXd fictitious_force_2(6);
        fictitious_force_2 << M_234 * skewSymmetric(Eigen::Vector3d::Zero()) * skewSymmetric(Eigen::Vector3d::Zero()) * (r_w_234 - r_w_2),
                              skewSymmetric(Eigen::Vector3d::Zero()) * (I_234 - M_234 * skewSymmetric(r_w_234 - r_w_2) * skewSymmetric(r_w_234 - r_w_2)) * Eigen::Vector3d::Zero();

        Eigen::VectorXd u_dot_7(6);
        u_dot_7 << Eigen::Vector3d::UnitX(),
                   Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_7 = spatial_inertia_2 * u_dot_7 + fictitious_force_2;

        Eigen::VectorXd s_7(6);
        s_7 << Eigen::Vector3d::Zero(),
               P_w_2;

        double mm7 = (s_7.transpose() * tau_7)(0);

        //8
        Eigen::VectorXd u_dot_8(6);
        u_dot_8 << Eigen::Vector3d::UnitY(),
                   Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_8 = spatial_inertia_2 * u_dot_8 + fictitious_force_2;

        Eigen::VectorXd s_8(6);
        s_8 << Eigen::Vector3d::Zero(),
               P_w_2;

        double mm8 = (s_8.transpose() * tau_8)(0);

        //9
        Eigen::VectorXd u_dot_9(6);
        u_dot_9 << Eigen::Vector3d::UnitZ(),
                   Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_9 = spatial_inertia_2 * u_dot_9 + fictitious_force_2;

        Eigen::VectorXd s_9(6);
        s_9 << Eigen::Vector3d::Zero(),
               P_w_2;

        double mm9 = (s_9.transpose() * tau_9)(0);

        //10
        Eigen::VectorXd u_dot_10(6);
        u_dot_10 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitX()) * (r_w_2 - r_w_b),
                    Eigen::Vector3d::UnitX();

        Eigen::VectorXd tau_10 = spatial_inertia_2 * u_dot_10 + fictitious_force_2;

        Eigen::VectorXd s_10(6);
        s_10 << Eigen::Vector3d::Zero(),
                P_w_2;

        double mm10 = (s_10.transpose() * tau_10)(0);

        //11
        Eigen::VectorXd u_dot_11(6);
        u_dot_11 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitY()) * (r_w_2 - r_w_b),
                    Eigen::Vector3d::UnitY();

        Eigen::VectorXd tau_11 = spatial_inertia_2 * u_dot_11 + fictitious_force_2;

        Eigen::VectorXd s_11(6);
        s_11 << Eigen::Vector3d::Zero(),
                P_w_2;

        double mm11 = (s_11.transpose() * tau_11)(0);

        //12
        Eigen::VectorXd u_dot_12(6);
        u_dot_12 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitZ()) * (r_w_2 - r_w_b),
                    Eigen::Vector3d::UnitZ();

        Eigen::VectorXd tau_12 = spatial_inertia_2 * u_dot_12 + fictitious_force_2;

        Eigen::VectorXd s_12(6);
        s_12 << Eigen::Vector3d::Zero(),
                P_w_2;

        double mm12 = (s_12.transpose() * tau_12)(0);

        //13
        Eigen::MatrixXd spatial_inertia_3(6, 6);
        spatial_inertia_3 << M_34 * Eigen::Matrix3d::Identity(), -M_34 * skewSymmetric(r_w_34 - r_w_3),
                             M_34 * skewSymmetric(r_w_34 - r_w_3), I_34 - M_34 * skewSymmetric(r_w_34 - r_w_3) * skewSymmetric(r_w_34 - r_w_3);

        Eigen::VectorXd fictitious_force_3(6);
        fictitious_force_3 << M_34 * skewSymmetric(Eigen::Vector3d::Zero()) * skewSymmetric(Eigen::Vector3d::Zero()) * (r_w_34 - r_w_3),
                              skewSymmetric(Eigen::Vector3d::Zero()) * (I_34 - M_34 * skewSymmetric(r_w_34 - r_w_3) * skewSymmetric(r_w_34 - r_w_3)) * Eigen::Vector3d::Zero();

        Eigen::VectorXd u_dot_13(6);
        u_dot_13 << Eigen::Vector3d::UnitX(),
                    Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_13 = spatial_inertia_3 * u_dot_13 + fictitious_force_3;

        Eigen::VectorXd s_13(6);
        s_13 << Eigen::Vector3d::Zero(),
                P_w_3;

        double mm13 = (s_13.transpose() * tau_13)(0);

        //14
        Eigen::VectorXd u_dot_14(6);
        u_dot_14 << Eigen::Vector3d::UnitY(),
                    Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_14 = spatial_inertia_3 * u_dot_14 + fictitious_force_3;

        Eigen::VectorXd s_14(6);
        s_14 << Eigen::Vector3d::Zero(),
                P_w_3;

        double mm14 = (s_14.transpose() * tau_14)(0);

        //15
        Eigen::VectorXd u_dot_15(6);
        u_dot_15 << Eigen::Vector3d::UnitZ(),
                    Eigen::Vector3d::Zero();

        Eigen::VectorXd tau_15 = spatial_inertia_3 * u_dot_15 + fictitious_force_3;

        Eigen::VectorXd s_15(6);
        s_15 << Eigen::Vector3d::Zero(),
                P_w_3;

        double mm15 = (s_15.transpose() * tau_15)(0);

        //16
        Eigen::VectorXd u_dot_16(6);
        u_dot_16 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitX()) * (r_w_3 - r_w_b),
                    Eigen::Vector3d::UnitX();

        Eigen::VectorXd tau_16 = spatial_inertia_3 * u_dot_16 + fictitious_force_3;

        Eigen::VectorXd s_16(6);
        s_16 << Eigen::Vector3d::Zero(),
                P_w_3;

        double mm16 = (s_16.transpose() * tau_16)(0);

        //17
        Eigen::VectorXd u_dot_17(6);
        u_dot_17 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitY()) * (r_w_3 - r_w_b),
                    Eigen::Vector3d::UnitY();

        Eigen::VectorXd tau_17 = spatial_inertia_3 * u_dot_17 + fictitious_force_3;

        Eigen::VectorXd s_17(6);
        s_17 << Eigen::Vector3d::Zero(),
                P_w_3;

        double mm17 = (s_17.transpose() * tau_17)(0);

        //18
        Eigen::VectorXd u_dot_18(6);
        u_dot_18 << Eigen::Vector3d::Zero() + skewSymmetric(Eigen::Vector3d::UnitZ()) * (r_w_3 - r_w_b),
                    Eigen::Vector3d::UnitZ();

        Eigen::VectorXd tau_18 = spatial_inertia_3 * u_dot_18 + fictitious_force_3;

        Eigen::VectorXd s_18(6);
        s_18 << Eigen::Vector3d::Zero(),
                P_w_3;

        double mm18 = (s_18.transpose() * tau_18)(0);

        mass_matrix_1 << mm1, mm2, mm3, mm4, mm5, mm6,
                         mm7, mm8, mm9, mm10, mm11, mm12,
                         mm13, mm14, mm15, mm16, mm17, mm18;
    }

    Eigen::MatrixXd mass_matrix_2(3, 3);
    {
        mass_matrix_2.setZero();

        //1
        Eigen::Vector3d r_w_12 = getCompPos(M_1, M_2, r_w_1COM, r_w_2COM);
        Eigen::Matrix3d I_12 = getCompInertia(M_1, M_2, I_1, I_2, r_w_1COM - r_w_12, r_w_2COM - r_w_12);
        double M_12 = M_1 + M_2;

        Eigen::Vector3d r_w_34 = getCompPos(M_3, M_4, r_w_3COM, r_w_4COM);
        Eigen::Matrix3d I_34 = getCompInertia(M_3, M_4, I_3, I_4, r_w_3COM - r_w_34, r_w_4COM - r_w_34);
        double M_34 = M_3 + M_4;

        Eigen::Vector3d r_w_1234 = getCompPos(M_12, M_34, r_w_12, r_w_34);
        Eigen::Matrix3d I_1234 = getCompInertia(M_12, M_34, I_12, I_34, r_w_12 - r_w_1234, r_w_34 - r_w_1234);
        double M_1234 = M_12 + M_34;

        Eigen::MatrixXd spatial_inertia_1(6, 6);
        spatial_inertia_1 << M_1234 * Eigen::Matrix3d::Identity(), -M_1234 * skewSymmetric(r_w_1234 - r_w_1),
                             M_1234 * skewSymmetric(r_w_1234 - r_w_1), I_1234 - M_1234 * skewSymmetric(r_w_1234 - r_w_1) * skewSymmetric(r_w_1234 - r_w_1);

        Eigen::VectorXd fictitious_force_1(6);
        fictitious_force_1 << M_1234 * skewSymmetric(Eigen::Vector3d::Zero()) * skewSymmetric(Eigen::Vector3d::Zero()) * (r_w_1234 - r_w_1),
                              skewSymmetric(Eigen::Vector3d::Zero()) * (I_1234 - M_1234 * skewSymmetric(r_w_1234 - r_w_1) * skewSymmetric(r_w_1234 - r_w_1)) * Eigen::Vector3d::Zero();

        Eigen::VectorXd u_dot_1(6);
        u_dot_1 << Eigen::Vector3d::Zero(),
                   P_w_1;

        Eigen::VectorXd tau_1 = spatial_inertia_1 * u_dot_1 + fictitious_force_1;
        double mm1 = (u_dot_1.transpose() * tau_1)(0);

        //2
        Eigen::Vector3d r_w_234 = getCompPos(M_2, M_34, r_w_2COM, r_w_34);
        Eigen::Matrix3d I_234 = getCompInertia(M_2, M_34, I_2, I_34, r_w_2COM - r_w_234, r_w_34 - r_w_234);
        double M_234 = M_2 + M_34;

        Eigen::MatrixXd spatial_inertia_2(6, 6);
        spatial_inertia_2 << M_234 * Eigen::Matrix3d::Identity(), -M_234 * skewSymmetric(r_w_234 - r_w_2),
                             M_234 * skewSymmetric(r_w_234 - r_w_2), I_234 - M_234 * skewSymmetric(r_w_234 - r_w_2) * skewSymmetric(r_w_234 - r_w_2);

        Eigen::VectorXd fictitious_force_2(6);
        fictitious_force_2 << M_234 * skewSymmetric(Eigen::Vector3d::Zero()) * skewSymmetric(Eigen::Vector3d::Zero()) * (r_w_234 - r_w_2),
                              skewSymmetric(Eigen::Vector3d::Zero()) * (I_234 - M_234 * skewSymmetric(r_w_234 - r_w_2) * skewSymmetric(r_w_234 - r_w_2)) * Eigen::Vector3d::Zero();

        Eigen::VectorXd u_dot_2(6);
        u_dot_2 << Eigen::Vector3d::Zero() + skewSymmetric(P_w_1) * (r_w_2 - r_w_1),
                   P_w_1;

        Eigen::VectorXd tau_2 = spatial_inertia_2 * u_dot_2 + fictitious_force_2;

        Eigen::VectorXd s_2(6);
        s_2 << Eigen::Vector3d::Zero(),
               P_w_2;

        double mm2 = (s_2.transpose() * tau_2)(0);

        //3
        Eigen::MatrixXd spatial_inertia_3(6, 6);
        spatial_inertia_3 << M_34 * Eigen::Matrix3d::Identity(), -M_34 * skewSymmetric(r_w_34 - r_w_3),
                             M_34 * skewSymmetric(r_w_34 - r_w_3), I_34 - M_34 * skewSymmetric(r_w_34 - r_w_3) * skewSymmetric(r_w_34 - r_w_3);

        Eigen::VectorXd fictitious_force_3(6);
        fictitious_force_3 << M_34 * skewSymmetric(Eigen::Vector3d::Zero()) * skewSymmetric(Eigen::Vector3d::Zero()) * (r_w_34 - r_w_3),
                              skewSymmetric(Eigen::Vector3d::Zero()) * (I_34 - M_34 * skewSymmetric(r_w_34 - r_w_3) * skewSymmetric(r_w_34 - r_w_3)) * Eigen::Vector3d::Zero();

        Eigen::VectorXd u_dot_3(6);
        u_dot_3 << Eigen::Vector3d::Zero() + skewSymmetric(P_w_1) * (r_w_3 - r_w_1),
                   P_w_1;

        Eigen::VectorXd tau_3 = spatial_inertia_3 * u_dot_3 + fictitious_force_3;

        Eigen::VectorXd s_3(6);
        s_3 << Eigen::Vector3d::Zero(),
               P_w_3;

        double mm3 = (s_3.transpose() * tau_3)(0);

        //4
        Eigen::MatrixXd spatial_inertia_4 = spatial_inertia_2;

        Eigen::MatrixXd fictitious_force_4 = fictitious_force_2;

        Eigen::VectorXd u_dot_4(6);
        u_dot_4 << Eigen::Vector3d::Zero(),
                   P_w_2;

        Eigen::VectorXd tau_4 = spatial_inertia_4 * u_dot_4 + fictitious_force_4;

        Eigen::VectorXd s_4(6);
        s_4 << Eigen::Vector3d::Zero(),
               P_w_2;

        double mm4 = (s_4.transpose() * tau_4)(0);

        //5
        Eigen::MatrixXd spatial_inertia_5 = spatial_inertia_3;

        Eigen::MatrixXd fictitious_force_5 = fictitious_force_3;

        Eigen::VectorXd u_dot_5(6);
        u_dot_5 << Eigen::Vector3d::Zero() + skewSymmetric(P_w_2) * (r_w_3 - r_w_2),
                   P_w_2;

        Eigen::VectorXd tau_5 = spatial_inertia_5 * u_dot_5 + fictitious_force_5;

        Eigen::VectorXd s_5(6);
        s_5 << Eigen::Vector3d::Zero(),
               P_w_2;

        double mm5 = (s_5.transpose() * tau_5)(0);

        //6
        Eigen::MatrixXd spatial_inertia_6 = spatial_inertia_3;

        Eigen::MatrixXd fictitious_force_6 = fictitious_force_3;

        Eigen::VectorXd u_dot_6(6);
        u_dot_6 << Eigen::Vector3d::Zero(),
                   P_w_3;

        Eigen::VectorXd tau_6 = spatial_inertia_6 * u_dot_6 + fictitious_force_6;

        double mm6 = (u_dot_6.transpose() * tau_6)(0);

        mass_matrix_2 << mm1, mm2, mm3,
                         mm2, mm4, mm5,
                         mm3, mm5, mm6;
    }

    Eigen::MatrixXd mass_matrix(9, 9);
    mass_matrix << spatial_inertia_matrix, mass_matrix_1.transpose(),
                   mass_matrix_1,          mass_matrix_2;
    
    return mass_matrix;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearitiesUsingRNE (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

    Eigen::Matrix3d R_w_b = toRot(gc[3], gc[4], gc[5], gc[6]);
    Eigen::Matrix3d R_b_1 = toRot(anymal.joint_hip.rpy);
    Eigen::Matrix3d R_1_2 = toRot(anymal.joint_thigh.rpy);
    Eigen::Matrix3d R_2_3 = toRot(anymal.joint_shank.rpy);
    Eigen::Matrix3d R_3_4 = toRot(anymal.joint_adapter.rpy);

    Eigen::Matrix3d R_1_1 = toRot(gc[7], Rotation::roll);
    Eigen::Matrix3d R_2_2 = toRot(gc[8], Rotation::pitch);
    Eigen::Matrix3d R_3_3 = toRot(gc[9], Rotation::pitch);

    Eigen::Matrix3d R_w_1 = R_w_b * R_b_1 * R_1_1;
    Eigen::Matrix3d R_w_2 = R_w_1 * R_1_2 * R_2_2;
    Eigen::Matrix3d R_w_3 = R_w_2 * R_2_3 * R_3_3;
    Eigen::Matrix3d R_w_4 = R_w_3 * R_3_4;

    Eigen::Vector3d r_b_bCOM = R_w_b * anymal.link_base.xyz;
    Eigen::Vector3d r_1_1COM = R_w_1 * anymal.link_hip.xyz;
    Eigen::Vector3d r_2_2COM = R_w_2 * anymal.link_thigh.xyz;
    Eigen::Vector3d r_3_3COM = R_w_3 * anymal.link_shank.xyz;
    Eigen::Vector3d r_4_4COM = R_w_4 * anymal.link_adapter.xyz;

    Eigen::Vector3d r_w_b = {gc[0], gc[1], gc[2]};
    Eigen::Vector3d r_b_1 = R_w_b * anymal.joint_hip.xyz;
    Eigen::Vector3d r_1_2 = R_w_1 * anymal.joint_thigh.xyz;
    Eigen::Vector3d r_2_3 = R_w_2 * anymal.joint_shank.xyz;
    Eigen::Vector3d r_3_4 = R_w_3 * anymal.joint_adapter.xyz;
    Eigen::Vector3d r_4_e = R_w_4 * anymal.joint_foot.xyz;

    Eigen::Vector3d r_w_bCOM = r_w_b + r_b_bCOM;
    Eigen::Vector3d r_b_1COM = r_b_1 + r_1_1COM;
    Eigen::Vector3d r_1_2COM = r_1_2 + r_2_2COM;
    Eigen::Vector3d r_2_3COM = r_2_3 + r_3_3COM;
    Eigen::Vector3d r_3_4COM = r_3_4 + r_4_4COM;

    Eigen::Vector3d r_w_1 = r_w_b + r_b_1;
    Eigen::Vector3d r_w_2 = r_w_1 + r_1_2;
    Eigen::Vector3d r_w_3 = r_w_2 + r_2_3;
    Eigen::Vector3d r_w_4 = r_w_3 + r_3_4;
    Eigen::Vector3d r_w_e = r_w_4 + r_4_e;

    Eigen::Vector3d r_w_1COM = r_w_b + r_b_1COM;
    Eigen::Vector3d r_w_2COM = r_w_1 + r_1_2COM;
    Eigen::Vector3d r_w_3COM = r_w_2 + r_2_3COM;
    Eigen::Vector3d r_w_4COM = r_w_3 + r_3_4COM;

    Eigen::Vector3d r_3_e = r_3_4 + r_4_e;
    Eigen::Vector3d r_2_e = r_2_3 + r_3_e;
    Eigen::Vector3d r_1_e = r_1_2 + r_2_e;
    Eigen::Vector3d r_b_e = r_b_1 + r_1_e;

    Eigen::Vector3d P_w_b = R_w_b * anymal.joint_base.axis;
    Eigen::Vector3d P_w_1 = R_w_1 * anymal.joint_hip.axis;
    Eigen::Vector3d P_w_2 = R_w_2 * anymal.joint_thigh.axis;
    Eigen::Vector3d P_w_3 = R_w_3 * anymal.joint_shank.axis;

    Eigen::Matrix3d I_b = R_w_b * toInertiaMatrix(anymal.link_base.inertia) * R_w_b.transpose();
    Eigen::Matrix3d I_1 = R_w_1 * toInertiaMatrix(anymal.link_hip.inertia) * R_w_1.transpose();
    Eigen::Matrix3d I_2 = R_w_2 * toInertiaMatrix(anymal.link_thigh.inertia) * R_w_2.transpose();
    Eigen::Matrix3d I_3 = R_w_3 * toInertiaMatrix(anymal.link_shank.inertia) * R_w_3.transpose();
    Eigen::Matrix3d I_4 = R_w_4 * toInertiaMatrix(anymal.link_adapter.inertia) * R_w_4.transpose();

    double M_b = anymal.link_base.mass;
    double M_1 = anymal.link_hip.mass;
    double M_2 = anymal.link_thigh.mass;
    double M_3 = anymal.link_shank.mass;
    double M_4 = anymal.link_adapter.mass;

    Eigen::Vector3d ang_vel_b = {gv[3], gv[4], gv[5]};
    Eigen::Vector3d ang_vel_1 = ang_vel_b + P_w_1 * gv[6];
    Eigen::Vector3d ang_vel_2 = ang_vel_1 + P_w_2 * gv[7];
    Eigen::Vector3d ang_vel_3 = ang_vel_2 + P_w_3 * gv[8];

    Eigen::Vector3d pos_vel_b = {gv[0], gv[1], gv[2]};
    Eigen::Vector3d pos_vel_1 = pos_vel_b + skewSymmetric(ang_vel_b) * r_b_1;
    Eigen::Vector3d pos_vel_2 = pos_vel_1 + skewSymmetric(ang_vel_1) * r_1_2;
    Eigen::Vector3d pos_vel_3 = pos_vel_2 + skewSymmetric(ang_vel_2) * r_2_3;

    Eigen::Vector3d gravity = {0, 0, 9.81};
    Eigen::VectorXd ga(9);
    ga << gravity,
          Eigen::Vector3d::Zero(),
          Eigen::Vector3d::Zero();

    Eigen::Vector3d ang_acc_b = {ga[3], ga[4], ga[5]};
    Eigen::Vector3d ang_acc_1 = ang_acc_b + skewSymmetric(ang_vel_b) * P_w_1 * gv[6] + P_w_1 * ga[6];
    Eigen::Vector3d ang_acc_2 = ang_acc_1 + skewSymmetric(ang_vel_1) * P_w_2 * gv[7] + P_w_2 * ga[7];
    Eigen::Vector3d ang_acc_3 = ang_acc_2 + skewSymmetric(ang_vel_2) * P_w_3 * gv[8] + P_w_3 * ga[8];

    Eigen::Vector3d pos_acc_b = {ga[0], ga[1], ga[2]};
    Eigen::Vector3d pos_acc_1 = pos_acc_b + skewSymmetric(ang_acc_b) * r_b_1 + skewSymmetric(ang_vel_b) * skewSymmetric(ang_vel_b) * r_b_1;
    Eigen::Vector3d pos_acc_2 = pos_acc_1 + skewSymmetric(ang_acc_1) * r_1_2 + skewSymmetric(ang_vel_1) * skewSymmetric(ang_vel_1) * r_1_2;
    Eigen::Vector3d pos_acc_3 = pos_acc_2 + skewSymmetric(ang_acc_2) * r_2_3 + skewSymmetric(ang_vel_2) * skewSymmetric(ang_vel_2) * r_2_3;

    //9
    Eigen::Vector3d r_w_34 = getCompPos(M_3, M_4, r_w_3COM, r_w_4COM);
    Eigen::Matrix3d I_34 = getCompInertia(M_3, M_4, I_3, I_4, r_w_3COM - r_w_34, r_w_4COM - r_w_34);
    double M_34 = M_3 + M_4;

    Eigen::MatrixXd spatial_inertia_9(6, 6);
    spatial_inertia_9 << M_34 * Eigen::Matrix3d::Identity(), -M_34 * skewSymmetric(r_w_34 - r_w_3),
                         M_34 * skewSymmetric(r_w_34 - r_w_3), I_34 - M_34 * skewSymmetric(r_w_34 - r_w_3) * skewSymmetric(r_w_34 - r_w_3);

    Eigen::VectorXd fictitious_force_9(6);
    fictitious_force_9 << M_34 * skewSymmetric(ang_vel_3) * skewSymmetric(ang_vel_3) * (r_w_34 - r_w_3),
                          skewSymmetric(ang_vel_3) * (I_34 - M_34 * skewSymmetric(r_w_34 - r_w_3) * skewSymmetric(r_w_34 - r_w_3)) * ang_vel_3;

    Eigen::VectorXd u_dot_9(6);
    u_dot_9 << pos_acc_3,
               ang_acc_3;

    Eigen::VectorXd tau_9 = spatial_inertia_9 * u_dot_9 + fictitious_force_9;

    Eigen::VectorXd s_9(6);
    s_9 << Eigen::Vector3d::Zero(),
           P_w_3;

    double b9 = (s_9.transpose() * tau_9)(0);

    //8
    Eigen::MatrixXd spatial_inertia_8(6, 6);
    spatial_inertia_8 << M_2 * Eigen::Matrix3d::Identity(), -M_2 * skewSymmetric(r_w_2COM - r_w_2),
                         M_2 * skewSymmetric(r_w_2COM - r_w_2), I_2 - M_2 * skewSymmetric(r_w_2COM - r_w_2) * skewSymmetric(r_w_2COM - r_w_2);

    Eigen::VectorXd fictitious_force_8(6);
    fictitious_force_8 << M_2 * skewSymmetric(ang_vel_2) * skewSymmetric(ang_vel_2) * (r_w_2COM - r_w_2),
                          skewSymmetric(ang_vel_2) * (I_2 - M_2 * skewSymmetric(r_w_2COM - r_w_2) * skewSymmetric(r_w_2COM - r_w_2)) * ang_vel_2;

    Eigen::VectorXd u_dot_8(6);
    u_dot_8 << pos_acc_2,
               ang_acc_2;

    Eigen::VectorXd reaction_force_9(6);
    reaction_force_9 << tau_9.head(3),
                        tau_9.tail(3) + skewSymmetric(r_2_3) * tau_9.head(3);

    Eigen::VectorXd tau_8 = spatial_inertia_8 * u_dot_8 + fictitious_force_8 + reaction_force_9;

    Eigen::VectorXd s_8(6);
    s_8 << Eigen::Vector3d::Zero(),
           P_w_2;

    double b8 = (s_8.transpose() * tau_8)(0);

    //7
    Eigen::MatrixXd spatial_inertia_7(6, 6);
    spatial_inertia_7 << M_1 * Eigen::Matrix3d::Identity(), -M_1 * skewSymmetric(r_w_1COM - r_w_1),
                         M_1 * skewSymmetric(r_w_1COM - r_w_1), I_1 - M_1 * skewSymmetric(r_w_1COM - r_w_1) * skewSymmetric(r_w_1COM - r_w_1);

    Eigen::VectorXd fictitious_force_7(6);
    fictitious_force_7 << M_1 * skewSymmetric(ang_vel_1) * skewSymmetric(ang_vel_1) * (r_w_1COM - r_w_1),
                          skewSymmetric(ang_vel_1) * (I_1 - M_1 * skewSymmetric(r_w_1COM - r_w_1) * skewSymmetric(r_w_1COM - r_w_1)) * ang_vel_1;

    Eigen::VectorXd u_dot_7(6);
    u_dot_7 << pos_acc_1,
               ang_acc_1;

    Eigen::VectorXd reaction_force_8(6);
    reaction_force_8 << tau_8.head(3),
                        tau_8.tail(3) + skewSymmetric(r_1_2) * tau_8.head(3);

    Eigen::VectorXd tau_7 = spatial_inertia_7 * u_dot_7 + fictitious_force_7 + reaction_force_8;

    Eigen::VectorXd s_7(6);
    s_7 << Eigen::Vector3d::Zero(),
           P_w_1;

    double b7 = (s_7.transpose() * tau_7)(0);

    //6
    Eigen::MatrixXd spatial_inertia_6(6, 6);
    spatial_inertia_6 << M_b * Eigen::Matrix3d::Identity(), -M_b * skewSymmetric(r_w_bCOM - r_w_b),
                         M_b * skewSymmetric(r_w_bCOM - r_w_b), I_b - M_b * skewSymmetric(r_w_bCOM - r_w_b) * skewSymmetric(r_w_bCOM - r_w_b);

    Eigen::VectorXd fictitious_force_6(6);
    fictitious_force_6 << M_b * skewSymmetric(ang_vel_b) * skewSymmetric(ang_vel_b) * (r_w_bCOM - r_w_b),
                          skewSymmetric(ang_vel_b) * (I_b - M_b * skewSymmetric(r_w_bCOM - r_w_b) * skewSymmetric(r_w_bCOM - r_w_b)) * ang_vel_b;

    Eigen::VectorXd u_dot_6(6);
    u_dot_6 << pos_acc_b,
               ang_acc_b;

    Eigen::VectorXd reaction_force_7(6);
    reaction_force_7 << tau_7.head(3),
                        tau_7.tail(3) + skewSymmetric(r_b_1) * tau_7.head(3);

    Eigen::VectorXd tau_6 = spatial_inertia_6 * u_dot_6 + fictitious_force_6 + reaction_force_7;

    Eigen::VectorXd s_6(6);
    s_6 << Eigen::Vector3d::Zero(),
           R_w_b * Eigen::Vector3d::UnitZ();

    double b6 = (s_6.transpose() * tau_6)(0);

    //5
    Eigen::VectorXd s_5(6);
    s_5 << Eigen::Vector3d::Zero(),
           R_w_b * Eigen::Vector3d::UnitY();

    double b5 = (s_5.transpose() * tau_6)(0);

    //4
    Eigen::VectorXd s_4(6);
    s_4 << Eigen::Vector3d::Zero(),
           R_w_b * Eigen::Vector3d::UnitX();

    double b4 = (s_4.transpose() * tau_6)(0);

    //3
    Eigen::VectorXd s_3(6);
    s_3 << R_w_b * Eigen::Vector3d::UnitZ(),
           Eigen::Vector3d::Zero();

    double b3 = (s_3.transpose() * tau_6)(0);

    //2
    Eigen::VectorXd s_2(6);
    s_2 << R_w_b * Eigen::Vector3d::UnitY(),
           Eigen::Vector3d::Zero();

    double b2 = (s_2.transpose() * tau_6)(0);

    //1
    Eigen::VectorXd s_1(6);
    s_1 << R_w_b * Eigen::Vector3d::UnitX(),
           Eigen::Vector3d::Zero();

    double b1 = (s_1.transpose() * tau_6)(0);

    Eigen::VectorXd nonlinearities(9);
    nonlinearities << b1, b2, b3, b4, b5, b6, b7, b8, b9;

    return nonlinearities;
}

