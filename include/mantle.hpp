#ifndef MANTLE_H
#define MANTLE_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <numeric>
#include <any>
#ifdef OPENMP
#include <omp.h>
#endif
#include "mantle_config.hpp"
#include "mantle_particle.hpp"
#include "mantle_particle_restricted.hpp"

typedef std::vector<Particle> (*userParticles)();
typedef std::vector<std::vector<ParticleHistory>> m_history;

struct reb_simulation;
class Mantle {

	private:
	reb_simulation* r;
	std::vector<Particle> particles;

    public:
	std::optional<m_coll_context> collisionContext;
    double MU;
    double LU;
    double TU;
	bool has_additional_forces;
	bool has_ode;

	public:
	Mantle(const std::vector<std::unique_ptr<struct MConfig>>&, std::vector<Particle>);
	~Mantle();
	m_history evalSim(double, int);

	public:
	long getCollisionLog() const;

	private:
	void addHash();
	void addAdditionalForce();
	void addParticles();
	void addCollision(m_coll_context);
	void addOde();

};

void m_additionalForcesPar(struct reb_simulation* const);

void m_additionalForces(struct reb_simulation* const);

int m_collisionResolve(struct reb_simulation*, struct reb_collision);

void m_ode(struct reb_ode* const, double* const, const double* const, const double);

#endif // MANTLE_H

