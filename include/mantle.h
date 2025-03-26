#include <vector>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <numeric>
#include <any>
#ifdef OPENMP
#include <omp.h>
#endif
#include "mantle_particle.h"
#include "mantle_particle_restricted.h"

typedef std::vector<Particle> (*userParticles)();
typedef std::vector<std::vector<ParticleHistory>> m_history;

enum class gravity_solvers {BASIC, COMPENSATED, TREE, JACOBI, NONE};
enum class collision_solvers {DIRECT, LINE, TREE, LINETREE, NONE};
enum class integrators {IAS15, WHFAST, LEAPFROG, SEI, NONE};

struct m_config {

	gravity_solvers gravity;
	collision_solvers collision;
	integrators integrator;
	double dt;
	int active;

	double opening_angle2;
	double softening;

	double box_size;
	int x_boxes;
	int y_boxes;
	int z_boxes;

	int threads;

	m_config() {
		gravity = gravity_solvers::BASIC;
		collision = collision_solvers::NONE;
		integrator = integrators::IAS15;
		opening_angle2 = 1.5;
		softening = 0.2;
		box_size = 1; x_boxes = 1; y_boxes = 1; z_boxes = 1;
		threads = 1;
		dt = -1;
		active = -1;
	}

};

struct m_sim_params {

	double G;
	double MU;
	double LU;
	double TU;

	m_config config;

	m_sim_params() : config() {G = 1.0; MU = 1.0; LU = 1.0; TU = 1.0;}

};

struct reb_simulation;
class Mantle {

	private:
	reb_simulation* r;
	m_sim_params reb_sim_params;
	std::optional<userParticles> userParticlesFn;
	std::vector<Particle> particles;
	std::optional<m_coll_context> collisionContext;
	bool has_additional_forces;
	bool has_ode;

	public:
	Mantle(m_sim_params, userParticles);
	Mantle(m_sim_params, std::vector<Particle>);
	~Mantle();
	void setAdditionalForce();
	void setCollisionContext(m_coll_context);
	void setOde();
	void initMantle();
	m_history evalSim(double, int);

	public:
	long getCollisionLog() const;

	private:
	void updateParams();
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

