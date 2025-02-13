#ifndef PARTICLE_RES_H
#define PARTICLE_RES_H

#include <variant>
#include <typeinfo>
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
	bool stateUpdated;
	Meta* meta_props;
	struct reb_particle* reb_p;

	public:
	ParticleRestricted();
	ParticleRestricted(struct reb_particle& p, double);
	void setState(double x, double y, double z, double vx, double vy, double vz);

	public:
	const std::any& getParticleProperties() const;
	const m_coll_context& getCollisionContext() const;
	const std::optional<std::vector<m_force_func>>& getAdditionalForces() const;

	public:
	void updateRebParticle();
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
	ParticleLite(struct reb_particle& p, double);
	void setState(double x, double y, double z, double vx, double vy, double vz);

	public:
	template <typename T>
	bool isType() const {return typeid(T) == p->getParticleProperties().type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(p->getParticleProperties());} else {return std::nullopt;}}

};

#endif // PARTICLE_RES_H

