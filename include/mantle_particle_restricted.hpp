#ifndef PARTICLE_RES_H
#define PARTICLE_RES_H

#include <variant>
#include <typeinfo>
#include <optional>
#include <vector>
#include <any>
#include "mantle_particle.hpp"

struct reb_particle;
class ParticleRestricted {

	public:
	m_cart_state_vector state;
	Meta& meta_props;

	public:
    ParticleRestricted(m_cart_state_vector _state, Meta& m) : meta_props(m) {state = _state;}
	ParticleRestricted& operator=(const ParticleRestricted& other) {
        if (this != &other) {
            meta_props = other.meta_props;
        }
        return *this;
    }

	public:
    const std::any& getParticleProperties() const {return meta_props.particle_properties;}
    const m_coll_context& getCollisionContext() const {return meta_props.collision_context;}
    const std::optional<std::vector<m_force_func>>& getAdditionalForces() const {return meta_props.additional_forces;}
    m_cart_state_vector updateRebParticle() {return state;}
	const std::optional<std::reference_wrapper<std::vector<Ode>>> getOdes() {return meta_props.getOdes();}
    std::optional<std::vector<double>> getOdey(uint32_t id) const {meta_props.syncOde(); return meta_props.getOdey(id);}
	uint32_t getHash() const {return meta_props.getHash();}
    bool isRemovable() const {return meta_props.isRemovable();}

	public:
	void incCollision() {meta_props.incCollision();}
	void dim();
	void nonDim();

};

class ParticleLite {

	private:
	ParticleRestricted& p;

	public:
	m_cart_state_vector& state;

	public:
    ParticleLite(ParticleRestricted &p_) : p(p_), state(p_.state) {}
	void setState(double x, double y, double z, double vx, double vy, double vz);

	public:
	template <typename T>
	bool isType() const {return typeid(T) == p.getParticleProperties().type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(p.getParticleProperties());} else {return std::nullopt;}}
    std::optional<std::vector<double>> getOdey(uint32_t id) const {return p.getOdey(id);}

};

class ParticleHistory {

	private:
	Meta meta_props;

	public:
	m_cart_state_vector state;

	public:
	ParticleHistory() {};
	ParticleHistory(m_cart_state_vector, Meta);

    public:
	bool collectHistory() const {return meta_props.collectHistory();}

	public:
	template <typename T>
	bool isType() const {return typeid(T) == meta_props.particle_properties.type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(meta_props.particle_properties);} else {return std::nullopt;}}
	int getCollisionLog() const {return meta_props.getCollisionLog();}
    std::optional<std::vector<double>> const getOdey(uint32_t id) {return meta_props.getOdey(id);};

	public:
	void dim();
	void nonDim();

};

#endif // PARTICLE_RES_H

