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
    ParticleRestricted(m_cart_state_vector _state, Meta& m) : meta_props(m) {state = _state; size = 0;}
	ParticleRestricted& operator=(const ParticleRestricted& other) {
        if (this != &other) {
            meta_props = other.meta_props;
        }
        return *this;
    }
    ~ParticleRestricted() {
        if (size > 0) {
            for (uint32_t i=0; i<size; i++) {
                Ode::nonDim(ys[i], 1.0, meta_props.TU);
            }
        }
    }

    private:
    uint32_t size;
    std::array<double*, 10> ys;

	public:
    const std::any& getParticleProperties() const {return meta_props.particle_properties;}
    const m_coll_context& getCollisionContext() const {return meta_props.collision_context;}
    const std::optional<std::vector<m_force_func>>& getAdditionalForces() const {return meta_props.additional_forces;}
    m_cart_state_vector updateRebParticle() {return state;}
	const std::optional<std::reference_wrapper<std::vector<Ode>>> getOdes() {return meta_props.getOdes();}
    std::optional<double*> getOdey(uint32_t id) {auto y = meta_props.getOdey(id); if (y.has_value()) {ys[size] = y.value(); size++;}; return y;}
	uint32_t getHash() const {return meta_props.getHash();}
    bool isRemovable() const {return meta_props.isRemovable();}

	public:
	void incCollision() {meta_props.incCollision();}
    inline void nonDim() {Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);dim_(state);}
    inline void dim() {Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU); dim_(state);}

};

class ParticleLite {

	private:
	ParticleRestricted& p;

	public:
	m_cart_state_vector& state;

	public:
    ParticleLite(ParticleRestricted &p_) : p(p_), state(p_.state) {}
    void setState(double x, double y, double z, double vx, double vy, double vz) {state.x = x; state.y = y; state.z = z; state.vx = vx; state.vy = vy; state.vz = vz;}

	public:
	template <typename T>
	bool isType() const {return typeid(T) == p.getParticleProperties().type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(p.getParticleProperties());} else {return std::nullopt;}}
    std::optional<double*> getOdey(uint32_t id) const {return p.getOdey(id);}

};

class ParticleHistory {

	private:
	Meta meta_props;

	public:
	m_cart_state_vector state;

	public:
	ParticleHistory() {};
    ParticleHistory(m_cart_state_vector _state, Meta m) {state = _state; meta_props = m; meta_props.syncOde();}

    public:
	bool collectHistory() const {return meta_props.collectHistory();}

	public:
	template <typename T>
	bool isType() const {return typeid(T) == meta_props.particle_properties.type();}
	template <typename T>
	std::optional<T> getParticleProperties() const {if (isType<T>()) {return std::any_cast<T>(meta_props.particle_properties);} else {return std::nullopt;}}
	int getCollisionLog() const {return meta_props.getCollisionLog();}
    std::optional<std::array<double, MAX_STATE>> const getOdey(uint32_t id) {return meta_props.getSaveOdey(id);};

	public:
    inline void nonDim() {Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU); dim_(state);}
    inline void dim() {Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU); dim_(state);}

};

#endif // PARTICLE_RES_H

