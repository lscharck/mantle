#include <vector>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <numeric>
#include <any>
#include "mantle_particle.h"
#include "mantle_particle_restricted.h"

typedef std::vector<Particle> (*userParticles)();

struct m_sim_params {

	double G;
	double MU;
	double LU;
	double TU;
	double box_size;

	bool gravity;
	bool collision;
	bool default_collision;

	bool tree_gravity;
	bool tree_collision;

};

struct reb_simulation;
class Mantle {

	private:
	reb_simulation* r;
	m_sim_params reb_sim_params;
	userParticles userParticlesFn;
	std::vector<Particle> particles;
	std::optional<m_coll_context> collisionContext;
	std::vector<std::pair<uint32_t, uint32_t>> ode_to_particle;
	bool has_additional_forces;
	bool has_ode;

	public:
	Mantle(m_sim_params, userParticles);
	void setAdditionalForce(void);
	void setCollisionContext(m_coll_context);
	void setOde(void);
	void initMantle(void);
	std::vector<std::vector<ParticleHistory>> evalSim(double, int);

	private:
	void updateParams(void);
	void addAdditionalForce(void);
	void addParticles(void);
	void addCollision(m_coll_context);
	void addOde(void);

};

void m_additionalForces(struct reb_simulation* const);

int m_collisionResolve(struct reb_simulation*, struct reb_collision);

void m_ode(struct reb_ode* const, double* const, const double* const, const double);

inline void dimState(m_cart_state_vector& p, double LU, double TU) {

    double VU = LU/TU;
    p.x = p.x * LU;
    p.y = p.y * LU;
    p.z = p.z * LU;
    p.vx = p.vx * VU;
    p.vy = p.vy * VU;
    p.vz = p.vz * VU;
    p.t = p.t * TU;

}

inline void nonDimState(m_cart_state_vector& p, double LU, double TU) {

    double VU = LU/TU;
    p.x = p.x / LU;
    p.y = p.y / LU;
    p.z = p.z / LU;
    p.vx = p.vx / VU;
    p.vy = p.vy / VU;
    p.vz = p.vz / VU;
    p.t = p.t / TU;

}

