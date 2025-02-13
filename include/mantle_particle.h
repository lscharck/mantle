#ifndef PARTICLE_H
#define PARTICLE_H

#include <variant>
#include <cstdint>
#include <optional>
#include <vector>
#include <functional>
#include <any>

struct m_accel {

	double ax;
	double ay;
	double az;

	friend m_accel operator+(const m_accel& a, const m_accel& b);
	friend m_accel operator*=(double s, const m_accel& a);
	friend m_accel operator*=(const m_accel& a, double s);

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
typedef void (*m_ode_func)(double* const, const ParticleLite&);

class Ode {

	public:
	m_ode_func ode_func;
	std::vector<double> y;
	std::any user_properties;
	uint32_t id;

};

class Meta {

	public:
	int collisions;
	bool count_collisions;
	bool collect_history;
	std::any user_properties;
	uint32_t hash;

	std::optional<std::vector<m_force_func>> additional_forces;
	m_coll_context collisionContext;
	std::vector<Ode> odes;

	double LU;
	double TU;
	double MU;

	public:
	Meta();
	Meta(const Meta* m);

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
	void setCollisionTracker();
	void setHistoryTracker();
	void setAdditionalForces(std::vector<m_force_func> function_array);
	void setODE(m_ode_func func, std::vector<double> y);

	template <typename T>
	void setParticleProperties(const T& properties) {meta_props.user_properties = properties;}

	protected:
	void setCollisionContext(m_coll_context collision_context);
	void setChar(double, double, double);
	void setHash(uint32_t);

	protected:
	const std::variant<m_cart_state_vector, m_kep_state_vector>& getState() const;
	const Meta& getMeta() const;
	const std::vector<Ode>& getODE() const;
	double getMass() const;
	double getRadius() const;
	uint32_t getHash() const;

	protected:
	void dim();
	void nonDim();

};

#endif // PARTICLE_H

