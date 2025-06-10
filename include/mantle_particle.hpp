#ifndef PARTICLE_H
#define PARTICLE_H

#include <variant>
#include <cstdint>
#include <optional>
#include <vector>
#include <functional>
#include <any>
#include <iostream>
#include <fstream>

#define MAX_STATE 10

struct m_accel {

	double ax;
	double ay;
	double az;

    friend m_accel operator+(const m_accel& t, const m_accel& u) {return {.ax = t.ax + u.ax, .ay = t.ay + u.ay, .az = t.az + u.az};}
    friend m_accel operator*(double s, const m_accel& u) {return {.ax = s*u.ax, .ay = s*u.ay, .az = s*u.az};}
    friend m_accel operator*(const m_accel& u, double s) {return {.ax = s*u.ax, .ay = s*u.ay, .az = s*u.az};}

};

struct m_cart_state_vector {

	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double t;

};

struct m_kep_state_vector {

	double a;
	double e;
	double inc;
	double Omega;
	double omega;
	double f;
	double t;

};

template <typename Op>
struct Scale {

    Op op;

    Scale(Op o, double l, double t) : op(o), LU(l), TU(t) {}

    void operator()(m_cart_state_vector& state) {

        state.x = op(state.x, LU);
        state.y = op(state.y, LU);
        state.z = op(state.z, LU);
        state.vx = op(state.vx, (LU / TU));
        state.vy = op(state.vy, (LU / TU));
        state.vz = op(state.vz, (LU / TU));
        state.t = op(state.t, TU);

    }

    void operator()(m_kep_state_vector& state) {

        state.a = op(state.a, LU);
        state.t = op(state.t, TU);

    }

    double LU;
    double TU;

};

class ParticleLite;
typedef m_accel (*m_force_func)(const ParticleLite&);
typedef int (*m_coll_func)(ParticleLite&, ParticleLite&);
typedef m_coll_func (*m_coll_context)(const ParticleLite&, const ParticleLite&);
typedef void (*m_ode_func)(double* const, const double* const, const ParticleLite&);

struct reb_ode;
class Ode {

	public:
	m_ode_func ode_func;
	std::array<double, MAX_STATE> init_y;
    uint32_t size;
	struct reb_ode* ho;
    std::array<double, MAX_STATE> save_y;
	std::any particle_properties;
	uint32_t id;
	uint32_t hash;
    double TU;

    public:
    void setOdey(struct reb_ode* _ho) {ho = _ho;};

    public:
    template <typename Op>
    inline static void scale(Op op, double* y, double s1, double s2) {

        for (uint32_t i=0; i<MAX_STATE; i+=2) {
            y[i] = op(y[i], s1);
        }

        for (uint32_t i=1; i<MAX_STATE; i+=2) {
            y[i] = op(y[i], s2);
        }

    }
    inline static void dim(double* y, double s1, double s2) {scale((std::divides<double>()), y, s1, s2);}
    inline static void nonDim(double* y, double s1, double s2) {scale((std::multiplies<double>()), y, s1, s2);}

};

class Meta {

	public:
	int collisions;
	bool count_collisions;
	bool collect_history;
    bool removable;
	std::any particle_properties;
	uint32_t hash;

	std::optional<std::vector<m_force_func>> additional_forces;
	m_coll_context collision_context;
	std::optional<std::vector<Ode>> odes;

	double LU;
	double TU;
	double MU;

	public:
	Meta();

    public:
   	void setOde(m_ode_func func, std::array<double, MAX_STATE> y, uint32_t id);
	void setCollisionTracker(void) {count_collisions = true;}
	void setHistoryTracker(void) {collect_history = true;}
	void setAdditionalForces(std::vector<m_force_func> function_array) {additional_forces = function_array;}
	void setCollisionContext(m_coll_context _collision_context) {collision_context = _collision_context;}
	void setChar(double _LU, double _TU, double _MU) {LU = _LU; TU = _TU; MU = _MU;}
	void setHash(uint32_t);
    void setIsRemovable(bool s) {removable = s;}
    void syncOde();

   	template <typename T>
	void setParticleProperties(const T& properties) {particle_properties = properties;}

	public:
	void incCollision() {if (count_collisions) {collisions++;}}
	int getCollisionLog() const {return collisions;}
	bool collectHistory() const {return collect_history;}
	const std::optional<std::reference_wrapper<std::vector<Ode>>> getOdes() {if (odes.has_value()) {return std::ref(odes.value());} else {return std::nullopt;}}
    const std::optional<double*> getOdey(uint32_t id) const;
    const std::optional<std::array<double, MAX_STATE>> getSaveOdey(uint32_t id) const;
	uint32_t getHash() const {return hash;}
    bool isRemovable() const {return removable;}

};

struct reb_particle;
class Mantle;
class Particle {

	friend class Mantle;

	private:
	std::variant<m_cart_state_vector, m_kep_state_vector> state;
	Meta meta_props;
	double m;
	double r;

	public:
	Particle(double mass, double radius);
	void setStateCart(double x, double y, double z, double vx, double vy, double vz);
	void setStateKep(double a, double e, double inc, double Omega, double omega, double f);
	void setCollisionTracker(void) {meta_props.setCollisionTracker();}
	void setHistoryTracker(void) {meta_props.setHistoryTracker();}
	void setAdditionalForces(std::vector<m_force_func> function_array) {meta_props.setAdditionalForces(function_array);}
	void setOde(m_ode_func func, std::array<double, MAX_STATE> y, uint32_t id) {meta_props.setOde(func, y, id);}
    void setIsRemovable(bool s) {meta_props.setIsRemovable(s);}

	template <typename T>
	void setParticleProperties(const T& properties) {meta_props.setParticleProperties(properties);}

	protected:
	void setCollisionContext(m_coll_context collision_context) {meta_props.setCollisionContext(collision_context);}
	void setChar(double LU, double TU, double MU) {meta_props.setChar(LU, TU, MU);}
	void setHash(uint32_t id) {meta_props.setHash(id);}

	protected:
	const std::variant<m_cart_state_vector, m_kep_state_vector>& getState() const {return state;}
	const Meta& getMeta() const {return meta_props;}
	double getMass() const {return m;}
	double getRadius() const {return r;}
	const std::optional<std::reference_wrapper<std::vector<Ode>>> getOdes() {return meta_props.getOdes();}
	uint32_t getHash() const {return meta_props.getHash();}

	protected:
    inline void nonDim() {
        Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);
        std::visit(dim_, state); m = m / meta_props.MU; r = r / meta_props.LU;
    }
    inline void dim() {
        Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU);
        std::visit(dim_, state); m = m * meta_props.MU; r = r * meta_props.LU;
    }

};

#endif // PARTICLE_H

