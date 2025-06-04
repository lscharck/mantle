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

struct m_accel {

	double ax;
	double ay;
	double az;

	friend m_accel operator+(const m_accel& a, const m_accel& b);
	friend m_accel operator*(double s, const m_accel& a);
	friend m_accel operator*(const m_accel& a, double s);

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

class ParticleLite;
typedef m_accel (*m_force_func)(const ParticleLite&);
typedef int (*m_coll_func)(ParticleLite&, ParticleLite&);
typedef m_coll_func (*m_coll_context)(const ParticleLite&, const ParticleLite&);
typedef std::vector<double> (*m_ode_func)(std::vector<double>, const ParticleLite&);

struct reb_ode;
class Ode {

	public:
	m_ode_func ode_func;
	std::vector<double> init_y;
    uint32_t size;
	struct reb_ode* ho;
    std::vector<double> save_y;
	std::any particle_properties;
	uint32_t id;
	uint32_t hash;
    double TU;

    public:
    void setOdey(struct reb_ode* _ho) {ho = _ho;};

    public:
    static std::vector<double> dim(std::vector<double>, double, double);
    static std::vector<double> nonDim(std::vector<double>, double, double);

};

class LUTManager {

	public:
	std::vector<std::pair<std::fstream, int>> LUTs;

	public:
	void addLUT(std::string);
	const std::fstream& getLUT(int) const;

};

class Meta {

	public:
	int collisions;
	bool count_collisions;
	bool collect_history;
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
   	void setOde(m_ode_func func, std::vector<double> y, uint32_t id);
	void setCollisionTracker(void) {count_collisions = true;}
	void setHistoryTracker(void) {collect_history = true;}
	void setAdditionalForces(std::vector<m_force_func> function_array) {additional_forces = function_array;}
	void setCollisionContext(m_coll_context _collision_context) {collision_context = _collision_context;}
	void setChar(double _LU, double _TU, double _MU) {LU = _LU; TU = _TU; MU = _MU;}
	void setHash(uint32_t);
    void saveOde();

   	template <typename T>
	void setParticleProperties(const T& properties) {particle_properties = properties;}

	public:
	void incCollision() {if (count_collisions) {collisions++;}}
	int getCollisionLog() const {return collisions;}
	bool collectHistory() const {return collect_history;}
	const std::optional<std::reference_wrapper<std::vector<Ode>>> getOdes() {if (odes.has_value()) {return std::ref(odes.value());} else {return std::nullopt;}}
    const std::optional<std::vector<double>> getOdey(uint32_t id) const;
	uint32_t getHash() const {return hash;}

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
	void setOde(m_ode_func func, std::vector<double> y, uint32_t id) {meta_props.setOde(func, y, id);}

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
	void dim();
	void nonDim();

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

#endif // PARTICLE_H

