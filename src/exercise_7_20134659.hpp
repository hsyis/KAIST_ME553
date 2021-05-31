#pragma once

#include <map>

namespace util {

constexpr double INF = std::numeric_limits<double>::infinity();

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

Eigen::Matrix3d quatToRot(const Eigen::Matrix<double, 4, 1>& quat) {
    return toRot(quat[0], quat[1], quat[2], quat[3]);
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

} // namespace util

template <typename Rep = double>
class Clock {
public:
    Clock()
        : time_(0)
        , step_(1) {
    }

    Rep getTimestep() const {
        return step_;
    }

    Clock& setTimestep(const Rep ts) {
        step_ = ts;
        return *this;
    }
    
    Clock& increase() {
        time_ += step_;
        return *this;
    }

    Clock& operator++() {
        this->increase();
        return *this;
    }

    Clock operator++(int) {
        Clock tmp = *this;
        this->increase();
        return tmp;
    }

    operator Rep() const {
        return time_;
    }

private:
    Rep time_;
    Rep step_;
};

enum class ObjectType { GROUND, SPHERE };

class SingleBodyObject {
public:
    explicit SingleBodyObject(ObjectType objectType)
        : objectType_(objectType)
        , position_(0, 0, 0)
        , linearVelocity_(0, 0, 0)
        , angularVelocity_(0, 0, 0) {
    }

    virtual ~SingleBodyObject() = default;

    virtual void integrate(const double dt, const Eigen::Vector3d& gravity) = 0;

    ObjectType getObjectType() const {
        return objectType_;
    }

    Eigen::Vector3d getPosition() const {
        return position_;
    }

    SingleBodyObject& setPosition(const Eigen::Vector3d& pos) {
        position_ = pos;
        return *this;
    }

    SingleBodyObject& setPosition(const double x, const double y, const double z) {
        return setPosition(Eigen::Vector3d{x, y, z});
    }

    Eigen::Vector3d getLinearVelocity() const {
        return linearVelocity_;
    }

    SingleBodyObject& setLinearVelocity(const Eigen::Vector3d& vel) {
        linearVelocity_ = vel;
        return *this;
    }

    SingleBodyObject& setLinearVelocity(const double x, const double y, const double z) {
        return setLinearVelocity(Eigen::Vector3d{x, y, z});
    }

    Eigen::Vector3d getAngularVelocity() const {
        return angularVelocity_;
    }

    SingleBodyObject& setAngularVelocity(const Eigen::Vector3d& vel) {
        angularVelocity_ = vel;
        return *this;
    }

    SingleBodyObject& setAngularVelocity(const double x, const double y, const double z) {
        return setAngularVelocity(Eigen::Vector3d{x, y, z});
    }

    Eigen::Matrix3d getBodyInertia() const {
        return bodyInertia_;
    }

    SingleBodyObject& addContactCandidate(std::shared_ptr<SingleBodyObject> object) {
        contactCandidates_.push_back(object);
        return *this;
    }

    auto contactCandidatesBegin() {
        return contactCandidates_.begin();
    }

    auto contactCandidatesEnd() {
        return contactCandidates_.end();
    }

protected:
    ObjectType objectType_;
    double mass_;
    Eigen::Vector3d position_;
    Eigen::Vector3d linearVelocity_;
    Eigen::Vector3d angularVelocity_;
    Eigen::Matrix3d rotation_;
    Eigen::Matrix3d bodyInertia_;
    Eigen::Matrix3d worldInertia_;
    std::vector<std::shared_ptr<SingleBodyObject>> contactCandidates_;
};

class Ground final : public SingleBodyObject {
public:
    Ground(const double height)
        : SingleBodyObject(ObjectType::GROUND)
        , height_(height) {
        mass_ = std::numeric_limits<decltype(mass_)>::infinity();
    }

    double getHeight() const {
        return height_;
    }

    Ground& setHeight(const double height) {
        height_ = height;
        return *this;
    }

    void integrate(const double dt, const Eigen::Vector3d& gravity) {
    }

private:
    double height_;
};

class Sphere final : public SingleBodyObject {
public:
    explicit Sphere(const double radius)
        : SingleBodyObject(ObjectType::SPHERE)
        , radius_(radius) {
        constexpr double density = 1000; // the density of water (kg/m^3)
        mass_ =  4.0 / 3.0 * M_PI * std::pow(radius_, 3) * density;
        bodyInertia_ = 2.0 / 5.0 * mass_ * std::pow(radius_, 2) * Eigen::Matrix3d::Identity();
    }

    double getRadius() const {
        return radius_;
    }

    Sphere& setRadius(const double radius) {
        radius_ = radius;
        return *this;
    }

    Eigen::MatrixXd getPositionalJacobian() {
        Eigen::MatrixXd jacobian(3, 6);
        jacobian << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero();
        return jacobian;
    }

    Eigen::MatrixXd getAngularJacobian() {
        Eigen::MatrixXd jacobian(3, 6);
        jacobian << Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Identity();
        return jacobian;
    }

    Eigen::MatrixXd getContactJacobian(const Eigen::Vector3d& contactPoint) {
        Eigen::MatrixXd jacobian(3, 6);
        Eigen::Vector3d r_com_contact = contactPoint - position_;
        jacobian << Eigen::Matrix3d::Identity(), -util::skewSymmetric(r_com_contact);
        return jacobian;
    }

    Eigen::MatrixXd getMassMatrix() {
        Eigen::MatrixXd posJ = getPositionalJacobian();
        Eigen::MatrixXd angJ = getAngularJacobian();
        return posJ.transpose() * mass_ * posJ + angJ.transpose() * bodyInertia_ * angJ;
    }

    Eigen::VectorXd getNonlinearity(const Eigen::Vector3d& gravity) {
        Eigen::MatrixXd posJ = getPositionalJacobian();
        return posJ.transpose() * mass_ * gravity;
    }

    void integrate(const double dt, const Eigen::Vector3d& gravity) {
        Eigen::Vector3d linearAcceleration = Eigen::Vector3d::Zero();
        Eigen::Vector3d angularAcceleration = Eigen::Vector3d::Zero();

        auto object = std::find_if(contactCandidates_.begin(), contactCandidates_.end(),
            [](const auto& obj) { return obj->getObjectType() == ObjectType::SPHERE; });

        auto sphere = reinterpret_cast<Sphere*>(object->get());

        if ((position_ - sphere->getPosition()).norm() < (radius_ + sphere->getRadius())) {
            // contact with another sphere
            double radius2 = sphere->getRadius();
            Eigen::Vector3d position2 = sphere->getPosition();

            Eigen::Vector3d contactPoint = (position_ * radius2 + position2 * radius_) / (radius_ + radius2);

            Eigen::MatrixXd J_c1 = getContactJacobian(contactPoint);
            Eigen::MatrixXd J_c2 = sphere->getContactJacobian(contactPoint);

            Eigen::MatrixXd J_a1 = getAngularJacobian();
            Eigen::MatrixXd J_a2 = sphere->getAngularJacobian();

            Eigen::Matrix3d M_app_inv = J_c1 * getMassMatrix().inverse() * J_c1.transpose() + J_c2 * sphere->getMassMatrix().inverse() * J_c2.transpose();
            Eigen::Vector3d bImp = -dt * J_c1 * getMassMatrix().inverse() * (getNonlinearity(gravity) + J_a1.transpose() * angularVelocity_.cross(bodyInertia_ * angularVelocity_))
                                   +dt * J_c2 * sphere->getMassMatrix().inverse() * (sphere->getNonlinearity(gravity) + J_a2.transpose() * sphere->angularVelocity_.cross(sphere->getBodyInertia() * sphere->getAngularVelocity()));

            Eigen::Vector3d lambda = Eigen::Vector3d::Zero();
            Eigen::Vector3d lambdaZ = Eigen::Vector3d::Zero();
            Eigen::Vector3d lambdaT = Eigen::Vector3d::Zero();
            Eigen::Vector3d lambdaBefore = Eigen::Vector3d::Constant(1);
            
            double alpha = 0.8;
            double C_z = alpha / M_app_inv(2, 2);
            double C_t = alpha / std::fmax(M_app_inv(0, 0), M_app_inv(1, 1));

            Eigen::Vector3d axisZ = position_ - position2;
            Eigen::Vector3d unitZ = axisZ / axisZ.norm();

            Eigen::Vector3d vImpBefore = (linearVelocity_ + angularVelocity_.cross(contactPoint - position_))
                                         -(sphere->getLinearVelocity() + sphere->getAngularVelocity().cross(contactPoint - position2) );

            // projected gauss-seidel method
            while ((lambda - lambdaBefore).norm() > 1e-10) {
                lambdaBefore = lambda;

                Eigen::Vector3d vImpAfter = vImpBefore + M_app_inv * lambda + bImp;

                // prox_z
                double valueZ = vImpAfter.transpose() * unitZ;
                double lambdaZ_val = ((lambdaZ - (valueZ * unitZ) * C_z).transpose() * unitZ);
                if (lambdaZ_val < 0) {
                    lambdaZ_val = 0;
                }
                lambdaZ = lambdaZ_val * unitZ;

                // prox_t
                double valueT = (vImpAfter - valueZ * unitZ).norm();

                Eigen::Vector3d tmp;
                if (valueT < 1e-10) {
                    tmp = Eigen::Vector3d::Zero();
                } else {
                    Eigen::Vector3d unitT = (vImpAfter - valueZ * unitZ) / valueT;
                    tmp = valueT * unitT;
                }

                double valueZ_2 = (lambdaT - tmp * C_t).transpose() * unitZ;

                double lambdaT_val = ((lambdaT - tmp * C_t) - valueZ_2 * unitZ).norm();
                if (lambdaT_val < 1e-10) {
                    lambdaT = Eigen::Vector3d::Zero();
                } else {
                    Eigen::Vector3d unitT = ((lambdaT - tmp * C_t) - valueZ_2 * unitZ) / lambdaT_val;
                    lambdaT = lambdaT_val * unitT;
                }

                if (lambdaT.transpose() * vImpBefore > 0) {
                    lambdaT.setZero();
                } else {
                    Eigen::Vector3d unitZ_2 = vImpBefore / vImpBefore.norm();
                    double valueZ_3 = vImpBefore.transpose() * unitZ_2;
                    lambdaT = valueZ_3 * unitZ_2;
                }

                // coulomb's friction model
                constexpr double coeff = 0.8;
                if (lambdaT.norm() > 1e-10 && lambdaT.norm() > lambdaZ.norm() * coeff) {
                    lambdaT = lambdaT / lambdaT.norm() * lambdaZ.norm() * coeff;
                }

                lambda = lambdaZ + lambdaT;
            }

            linearAcceleration = getPositionalJacobian() * getMassMatrix().inverse() * (J_c1.transpose() * lambda / dt - getNonlinearity(gravity));
            angularAcceleration = getAngularJacobian() * getMassMatrix().inverse() * (J_c1.transpose() * lambda / dt - getNonlinearity(gravity));
        } else if (position_(2) + 1e-5 < radius_) {
            // contact with ground
            Eigen::Vector3d contactPoint = {position_(0), position_(1), 0};

            Eigen::MatrixXd J_c = getContactJacobian(contactPoint);
            Eigen::MatrixXd J_a = getAngularJacobian();
            Eigen::Matrix3d M_app_inv = J_c * getMassMatrix().inverse() * J_c.transpose();
            Eigen::Vector3d bImp = -dt * J_c * getMassMatrix().inverse() * (getNonlinearity(gravity) + J_a.transpose() * angularVelocity_.cross(bodyInertia_ * angularVelocity_));

            Eigen::Vector3d lambda = Eigen::Vector3d::Zero();
            Eigen::Vector3d lambdaZ = Eigen::Vector3d::Zero();
            Eigen::Vector3d lambdaT = Eigen::Vector3d::Zero();
            Eigen::Vector3d lambdaBefore = Eigen::Vector3d::Constant(1);
            
            double alpha = 0.8;
            double C_z = alpha / M_app_inv(2, 2);
            double C_t = alpha / std::fmax(M_app_inv(0, 0), M_app_inv(1, 1));

            Eigen::Vector3d axisZ = {0, 0, 1};
            Eigen::Vector3d unitZ = axisZ / axisZ.norm();

            Eigen::Vector3d vImpBefore = linearVelocity_ + angularVelocity_.cross(contactPoint - position_);

            // projected gauss-seidel method
            while ((lambda - lambdaBefore).norm() > 1e-10) {
                lambdaBefore = lambda;

                Eigen::Vector3d vImpAfter = vImpBefore + M_app_inv * lambda + bImp;

                // prox_z
                double valueZ = vImpAfter.transpose() * unitZ;
                double lambdaZ_val = ((lambdaZ - (valueZ * unitZ) * C_z).transpose() * unitZ);
                if (lambdaZ_val < 0) {
                    lambdaZ_val = 0;
                }
                lambdaZ = lambdaZ_val * unitZ;

                // prox_t
                double valueT = (vImpAfter - valueZ * unitZ).norm();

                Eigen::Vector3d tmp;
                if (valueT < 1e-10) {
                    tmp = Eigen::Vector3d::Zero();
                } else {
                    Eigen::Vector3d unitT = (vImpAfter - valueZ * unitZ) / valueT;
                    tmp = valueT * unitT;
                }

                double valueZ_2 = (lambdaT - tmp * C_t).transpose() * unitZ;

                double lambdaT_val = ((lambdaT - tmp * C_t) - valueZ_2 * unitZ).norm();
                if (lambdaT_val < 1e-10) {
                    lambdaT = Eigen::Vector3d::Zero();
                } else {
                    Eigen::Vector3d unitT = ((lambdaT - tmp * C_t) - valueZ_2 * unitZ) / lambdaT_val;
                    lambdaT = lambdaT_val * unitT;
                }

                if (lambdaT.transpose() * vImpBefore > 0) {
                    lambdaT.setZero();
                } else {
                    Eigen::Vector3d unitZ_2 = vImpBefore / vImpBefore.norm();
                    double valueZ_3 = vImpBefore.transpose() * unitZ_2;
                    lambdaT = valueZ_3 * unitZ_2;
                }

                // coulomb's friction model
                constexpr double coeff = 0.8;
                if (lambdaT.norm() > 1e-10 && lambdaT.norm() > lambdaZ.norm() * coeff) {
                    lambdaT = lambdaT / lambdaT.norm() * lambdaZ.norm() * coeff;
                }

                lambda = lambdaZ + lambdaT;
            }

            linearAcceleration = getPositionalJacobian() * getMassMatrix().inverse() * (J_c.transpose() * lambda / dt - getNonlinearity(gravity));
            angularAcceleration = getAngularJacobian() * getMassMatrix().inverse() * (J_c.transpose() * lambda / dt - getNonlinearity(gravity));
        } else {
            // no contact
            linearAcceleration = gravity;
        }

        linearVelocity_ += linearAcceleration * dt;
        angularVelocity_ += angularAcceleration * dt;
        position_ += linearVelocity_ * dt;
    }

private:
    double radius_;
};

class SphereWorld {
public:
    SphereWorld() : gravity_{0, 0, -9.81} {
    }

    void setTimestep(const double ts) {
        clock_.setTimestep(ts);
    }

    void integrate() {
        clock_++;
        for (auto& p : objects_) {
            auto obj = p.second;
            obj->integrate(clock_.getTimestep(), gravity_);
        }
    }

    Sphere& getSphere(const std::string& name) {
        auto search = objects_.find(name);
        if (search != objects_.end()) {
            return *reinterpret_cast<Sphere*>(objects_[name].get());
        } else {
            return addSphere(name, 0);
        }
    }

    Sphere& addSphere(const std::string& name, const double radius) {
        auto sphere = std::make_shared<Sphere>(radius);
        updateContactCandidate(sphere);
        objects_[name] = sphere;
        return *reinterpret_cast<Sphere*>(objects_[name].get());
    }

    Ground& addGround(const double height) {
        const std::string name = "ground";
        auto ground = std::make_shared<Ground>(height);
        updateContactCandidate(ground);
        objects_[name] = ground;
        return *reinterpret_cast<Ground*>(objects_[name].get());
    }

private:
    void updateContactCandidate(std::shared_ptr<SingleBodyObject> object) {
        for (auto it = objects_.begin(); it != objects_.end(); it++) {
            it->second->addContactCandidate(object);
            object->addContactCandidate(it->second);
        }
    }
    
private:
    Clock<> clock_;
    Eigen::Vector3d gravity_;
    std::map<std::string, std::shared_ptr<SingleBodyObject>> objects_;
};

class SimulationClass {
public:
    SimulationClass() {
        world_.addGround(0.0);
        world_.setTimestep(0.001);

        world_.addSphere("sphere1", 0.5).setPosition(0.0, 0.0, 0.5);
        world_.addSphere("sphere2", 0.7).setPosition(0.1, 0.1, 3.0);
    }

    void integrate() {
        /// your code here
        world_.integrate();
    }

    void setPosition(raisim::Visuals * sphere1, raisim::Visuals * sphere2) {
        /// your code here
        auto pos1 = world_.getSphere("sphere1").getPosition();
        auto pos2 = world_.getSphere("sphere2").getPosition();

        sphere1->setPosition(pos1);
        sphere2->setPosition(pos2);
    }

    /// add state variables here
    SphereWorld world_;
};
