#ifndef PARTICLE_RES_H
#define PARTICLE_RES_H

#include <variant>
#include <optional>
#include <vector>
#include <any>
#include "mantle_particle.h"

struct reb_particle;
class ParticleRestricted {

	public:
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double t;
	double LU;
	double TU;
	double VU;
	std::optional<m_cart_state_vector> state;
	Meta* meta_props;
	struct reb_particle* reb_p;

	public:
	ParticleRestricted();
	ParticleRestricted(struct reb_particle& p, double);
	void setState(double x, double y, double z, double vx, double vy, double vz);

	public:
	std::any getParticleProperties() const;
	m_coll_context getCollisionContext();
	std::optional<std::vector<m_force_func>> getAdditionalForces();

	public:
	template <typename Op>
	void scale(Op);
	void dim();
	void nonDim();

};

class ParticleLite {

	private:
	ParticleRestricted* p;

	public:
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double t;

	public:
	ParticleLite(ParticleRestricted& p);
	void setState(double x, double y, double z, double vx, double vy, double vz);

	public:
	template <typename T>
	std::optional<T> getParticleProperties() const {return std::any_cast<T>(p->getParticleProperties());}

};

#endif // PARTICLE_RES_H

