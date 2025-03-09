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
	m_cart_state_vector state;
	bool stateUpdated;
	Meta& meta_props;
	struct reb_particle* reb_p;

	public:
	ParticleRestricted(struct reb_particle& p, double);

	public:
	const std::any& getParticleProperties() const;
	const m_coll_context& getCollisionContext() const;
	const std::optional<std::vector<m_force_func>>& getAdditionalForces() const;

	public:
	void incCollision() {meta_props.incCollision();}
	void updateRebParticle();
	void dim();
	void nonDim();

};

class ParticleLite {

	private:
	ParticleRestricted& p;

	public:
	m_cart_state_vector& state;

	public:
	ParticleLite(ParticleRestricted& p);
	ParticleLite(struct reb_particle& p, double);
	void setState(double x, double y, double z, double vx, double vy, double vz);

	public:
	template <typename T>
	bool isType() const {return typeid(T) == p.getParticleProperties().type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(p.getParticleProperties());} else {return std::nullopt;}}

};

class ParticleHistory {

	private:
	Meta meta_props;

	public:
	m_cart_state_vector state;

	public:
	ParticleHistory();
	ParticleHistory(struct reb_particle& p, double);

	public:
	template <typename T>
	bool isType() const {return typeid(T) == meta_props.particle_properties.type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(meta_props.particle_properties);} else {return std::nullopt;}}
	int getCollisionLog() const {return meta_props.getCollisionLog();}
	bool collectHistory() const {return meta_props.collectHistory();}

	public:
	void dim();
	void nonDim();

};

#endif // PARTICLE_RES_H

