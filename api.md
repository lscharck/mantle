# Mantle C++ API Reference

**Mantle** is an object-oriented C++ wrapper around the **REBOUND** C N-body astronomical simulation library (`rebound.h`). It provides type-safe configuration structures, automatic non-dimensionalization using characteristic units, per-particle metadata attachment, custom force and collision callbacks, and coupled Ordinary Differential Equation (ODE) integration.

## Core Simulation Engine

### `Mantle`
**Header:** `mantle.hpp` | **Source:** `mantle.cpp`

#### Purpose & Wrapped C Components
`Mantle` is the central simulation manager that orchestrates the lifecycle, configuration, particle initialization, and time integration of a REBOUND simulation. It wraps the core `reb_simulation` C struct and manages calls to:
- `reb_simulation_create` and `reb_simulation_free` for simulation memory management.
- `reb_simulation_add` and `reb_particle_from_orbit` for populating `reb_particle` structs (sorting particles so the most massive body is added first at index `0`).
- `reb_ode_create` and `reb_ode_free` for coupling custom `reb_ode` systems to individual particles.
- `reb_simulation_integrate` and `reb_mpi_finalize` (when compiled with `MPI`) for stepping the simulation through time and recording state snapshots.
- REBOUND function hooks (`r->additional_forces`, `r->collision_resolve`, `r->heartbeat`) and the `r->extras` pointer (which stores a back-pointer to the `Mantle` instance).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `collisionContext` | `std::optional<m_coll_context>` | Optional factory callback that returns a collision resolution function (`m_coll_func`) for a given pair of colliding particles. |
| `MU` | `double` | Characteristic mass unit used to non-dimensionalize particle masses before passing them to REBOUND. |
| `LU` | `double` | Characteristic length unit used to non-dimensionalize positions, radii, and orbital semi-major axes. |
| `TU` | `double` | Characteristic time unit used to non-dimensionalize velocities, accelerations, simulation durations, and ODE states. |
| `has_additional_forces` | `bool` | Flag indicating whether custom per-particle non-gravitational forces (`r->additional_forces`) are enabled. |
| `force_is_velocity_dep` | `bool` | Flag indicating whether custom forces depend on particle velocities (maps to `r->force_is_velocity_dependent`). |
| `has_ode` | `bool` | Flag indicating whether coupled `reb_ode` differential equations are attached to particles in the simulation. |
| `removeParticleTest` | `std::optional<m_remove_particle>` | Optional predicate evaluated at each heartbeat step to test whether a particle should be removed from the simulation. |
| `simTerminate` | `std::optional<m_sim_terminate>` | Optional predicate evaluated at each heartbeat step (receiving the active particle count `r->N`) to stop the simulation early via `reb_simulation_stop`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `Mantle` (Constructor) | `const std::vector<std::unique_ptr<struct MConfig>>& user_configs`<br>`std::vector<Particle> userParticles` | — | Allocates a new `reb_simulation` via `reb_simulation_create()`, links `r->extras` to `this`, applies each configuration object via `updateSim`, assigns sequential hashes to particles, registers C callbacks for collisions, additional forces, ODEs, and heartbeats, non-dimensionalizes all particles (placing the most massive particle first), and adds them to the REBOUND simulation. |
| `~Mantle` (Destructor) | None | — | Iterates through all particles to free any allocated `reb_ode` pointers via `reb_ode_free`, then frees the underlying `reb_simulation` struct via `reb_simulation_free`. |
| `evalSim` | `double tf`<br>`int nt`<br>`bool progress` | `m_history` | Integrates the simulation from $t = 0$ to final dimensional time `tf` across `nt` intervals (`nt + 1` total snapshots, including $t = 0$) using `reb_simulation_integrate`. At each step, converts internal `reb_particle` states back to dimensional units and records `ParticleHistory` snapshots for particles with history tracking enabled. Optionally prints completion percentage to `std::cout` when `progress` is `true`, halts early if `r->status == REB_STATUS_USER`, and calls `reb_mpi_finalize(r)` if `MPI` is enabled. |
| `evalSim` | `double tf`<br>`int nt` | `m_history` | Overload that calls `evalSim(tf, nt, false)` with console progress output disabled. |
| `getCollisionLog` | None | `long` | Returns the cumulative number of resolved collisions recorded by the underlying REBOUND simulation (`r->collisions_log_n`). |

## Particle Classes & Data Structures

### `Particle`
**Header:** `mantle_particle.hpp` | **Source:** `mantle_particle.cpp`

#### Purpose & Wrapped C Components
`Particle` is the primary user-facing class used to configure an individual body before passing it into the `Mantle` constructor. It encapsulates the particle's mass, physical radius, initial kinematic state (either Cartesian coordinates or Keplerian orbital elements), and a `Meta` instance. During `Mantle` initialization, each `Particle` is non-dimensionalized and converted into a REBOUND `reb_particle` struct via either direct Cartesian assignment (`reb_simulation_add`) or orbital conversion (`reb_particle_from_orbit`), with its `Meta` object linked to `reb_particle.ap`.

#### Public Data Members
`Particle` has no public data members. Its internal state (`state`, `meta_props`, `m`, `r`) is `private`, and its simulation-setup accessors and mutators (`getState`, `getMeta`, `getMass`, `getRadius`, `getOdes`, `getHash`, `setCollisionContext`, `setChar`, `setHash`, `nonDim`, `dim`) are `protected` and accessed by `Mantle` via `friend class Mantle`.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `Particle` (Constructor) | `double mass`<br>`double radius` | — | Initializes a particle with the specified physical `mass` (`m`) and `radius` (`r`) and default-constructs its `Meta` metadata container. |
| `setStateCart` | `double x`, `double y`, `double z`<br>`double vx`, `double vy`, `double vz` | `void` | Sets the particle's initial state to a 3D Cartesian position and velocity vector (`m_cart_state_vector`) with $t = 0$. |
| `setStateKep` | `double a`, `double e`, `double inc`<br>`double Omega`, `double omega`, `double f` | `void` | Sets the particle's initial state to classical Keplerian orbital elements (`m_kep_state_vector`) relative to the primary body (`r->particles[0]`), which `Mantle` converts via `reb_particle_from_orbit`. |
| `setCollisionTracker` | `void` | `void` | Enables per-particle collision counting in the underlying `Meta` instance (`count_collisions = true`). |
| `setHistoryTracker` | `void` | `void` | Marks this particle to be included in the `m_history` state snapshots returned by `Mantle::evalSim` (`collect_history = true`). |
| `setAdditionalForces` | `std::vector<m_force_func> function_array` | `void` | Attaches a list of custom acceleration callback functions (`m_force_func`) to be evaluated for this particle during REBOUND force calculations. |
| `setOde` | `m_ode_func func`<br>`std::vector<double> y`<br>`uint32_t id` | `void` | Registers a custom coupled ODE system (`Ode`) with derivative callback `func`, initial state vector `y`, and user identifier `id`. |
| `setIsRemovable` | `bool s` | `void` | Sets whether the particle is eligible to be removed from the simulation by the heartbeat `removeParticleTest` predicate. |
| `setParticleProperties<T>` | `const T& properties` | `void` | Template method that stores arbitrary user-defined particle properties of type `T` inside `meta_props.particle_properties` (`std::any`). |
| `setParticlePropertiesGen` | `const PropertiesGen properties` | `void` | Stores a predefined `PropertiesGen` struct (mass, radius, thermal, and optical properties) inside `meta_props.particle_properties_gen`. |

### `Meta`
**Header:** `mantle_particle.hpp` | **Source:** `mantle_particle.cpp`

#### Purpose & Wrapped C Components
`Meta` stores per-particle metadata, physical unit scaling factors, custom properties, callback lists, and coupled `Ode` objects. A pointer to each particle's `Meta` instance is cast to `void*` and stored in the `ap` (additional pointer) field of the corresponding C `reb_particle` struct (`p.ap`), allowing C callbacks (`m_additionalForces`, `m_collisionResolve`, `m_ode`, `m_heartbeat`) to recover the C++ metadata for any particle in `r->particles`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `collisions` | `int` | Running count of collisions involving this particle (incremented when `count_collisions` is `true`, initialized to `0`). |
| `count_collisions` | `bool` | Flag indicating whether collision events should increment `collisions` (initialized to `false`). |
| `collect_history` | `bool` | Flag indicating whether `Mantle::evalSim` should record `ParticleHistory` snapshots for this particle (initialized to `false`). |
| `removable` | `bool` | Flag indicating whether the particle can be removed via `reb_simulation_remove_particle_by_hash` during heartbeat checks (initialized to `true`). |
| `particle_properties` | `std::any` | Type-erased container holding custom user-defined particle property structs. |
| `particle_properties_gen` | `std::optional<PropertiesGen>` | Optional container holding standard physical and thermal properties (`PropertiesGen`). |
| `hash` | `uint32_t` | Unique identifier assigned to the particle and stored in `reb_particle.hash` for lookup via `reb_simulation_particle_by_hash` (initialized to `0`). |
| `additional_forces` | `std::optional<std::vector<m_force_func>>` | Optional list of custom non-gravitational force callbacks attached to the particle. |
| `collision_context` | `m_coll_context` | Collision context dispatcher function assigned when collision handling is configured. |
| `odes` | `std::optional<std::vector<Ode>>` | Optional vector of coupled `Ode` instances attached to this particle. |
| `LU` | `double` | Characteristic length unit for scaling this particle (initialized to `1.0`). |
| `TU` | `double` | Characteristic time unit for scaling this particle (initialized to `1.0`). |
| `MU` | `double` | Characteristic mass unit for scaling this particle (initialized to `1.0`). |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `Meta` (Constructor) | None | — | Initializes `collisions = 0`, `count_collisions = false`, `collect_history = false`, `removable = true`, `hash = 0`, and characteristic units `LU = 1.0`, `TU = 1.0`, `MU = 1.0`. |
| `setOde` | `m_ode_func func`<br>`std::vector<double> y`<br>`uint32_t id` | `void` | Creates an `Ode` instance with callback `func`, copies initial state `y` into `init_y` and `save_y`, copies current particle properties, and appends it to `odes`. |
| `setCollisionTracker` | `void` | `void` | Enables collision counting (`count_collisions = true`). |
| `setHistoryTracker` | `void` | `void` | Enables state history collection (`collect_history = true`). |
| `setAdditionalForces` | `std::vector<m_force_func> function_array` | `void` | Assigns the vector of custom force callbacks to `additional_forces`. |
| `setCollisionContext` | `m_coll_context _collision_context` | `void` | Assigns the collision context selector callback to `collision_context`. |
| `setChar` | `double _LU`<br>`double _TU`<br>`double _MU` | `void` | Sets the characteristic length (`LU`), time (`TU`), and mass (`MU`) scaling units. |
| `setHash` | `uint32_t _hash` | `void` | Sets `hash` for the particle and propagates the hash value to all attached `Ode` instances in `odes`. |
| `setIsRemovable` | `bool s` | `void` | Sets the `removable` flag to `s`. |
| `syncOde` | None | `void` | Dimensionalizes each attached `reb_ode` state array (`ode.ho->y`), copies the values into `ode.save_y` for history logging, and restores `ode.ho->y` to non-dimensional units. |
| `nonDimId` | `uint32_t id` | `void` | Finds the attached `Ode` matching `id` and converts its live REBOUND `ho->y` array back to non-dimensional units via `Ode::nonDim`. |
| `setParticleProperties<T>` | `const T& properties` | `void` | Template method that assigns a custom property object of type `T` to `particle_properties`. |
| `setParticlePropertiesGen` | `const PropertiesGen properties` | `void` | Assigns a `PropertiesGen` struct to `particle_properties_gen`. |
| `incCollision` | None | `void` | Increments `collisions` by `1` if `count_collisions` is `true`. |
| `getCollisionLog` | None | `int` | Returns the total number of collisions recorded for this particle (`collisions`). |
| `collectHistory` | None | `bool` | Returns `true` if history collection is enabled for this particle (`collect_history`). |
| `getOdes` | None | `const std::optional<std::reference_wrapper<std::vector<Ode>>>` | Returns an optional reference wrapper to the `odes` vector if present, or `std::nullopt`. |
| `getOdey` | `uint32_t id` | `const std::optional<double*>` | Finds the `Ode` with the specified `id`, dimensionalizes its live REBOUND `ho->y` array in-place, and returns the raw `double*` pointer (or `std::nullopt` if not found). |
| `getSaveOdey` | `uint32_t id` | `const std::optional<std::array<double, MAX_STATE>>` | Returns the saved dimensionalized ODE state array (`save_y`) for the `Ode` matching `id` (or `std::nullopt` if not found). |
| `getHash` | None | `uint32_t` | Returns the particle's unique `hash` identifier. |
| `isRemovable` | None | `bool` | Returns whether the particle is marked as removable (`removable`). |

### `Ode`
**Header:** `mantle_particle.hpp` | **Source:** `mantle_particle.cpp`

#### Purpose & Wrapped C Components
`Ode` wraps REBOUND's `reb_ode` struct (`ho`), which allows arbitrary systems of ordinary differential equations to be integrated simultaneously with the N-body equations of motion. In `Mantle::addOde`, each `Ode` creates a `reb_ode` via `reb_ode_create(r, ode.size)`, binds `ho->derivatives = m_ode`, copies the non-dimensionalized initial state into `ho->y`, and sets `ho->ref` to point back to the `Ode` instance.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `ode_func` | `m_ode_func` | User-supplied callback function that computes the state derivatives `yDot` given state `y`, time `t`, and a `ParticleLite` view. |
| `init_y` | `std::array<double, MAX_STATE>` | Initial state vector (up to `MAX_STATE = 10` elements). |
| `size` | `uint32_t` | Number of active state variables in the ODE system. |
| `ho` | `struct reb_ode*` | Pointer to the underlying REBOUND `reb_ode` C structure. |
| `save_y` | `std::array<double, MAX_STATE>` | Cached dimensional state array synchronized during `ParticleHistory` snapshot creation. |
| `particle_properties` | `std::any` | Copy of the parent particle's custom properties (`std::any`) at the time the ODE was registered. |
| `particle_properties_gen` | `std::optional<PropertiesGen>` | Copy of the parent particle's `PropertiesGen` struct at the time the ODE was registered. |
| `id` | `uint32_t` | User-assigned identifier used to query this specific ODE system via `getOdey(id)`. |
| `hash` | `uint32_t` | Hash of the parent particle, used by `m_ode` to locate the `reb_particle` via `reb_simulation_particle_by_hash`. |
| `TU` | `double` | Characteristic time unit used to scale even-indexed and odd-indexed state and derivative components. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `setOdey` | `struct reb_ode* _ho` | `void` | Stores the pointer to the allocated REBOUND `reb_ode` struct in `ho`. |
| `scale<Op>` (Static) | `Op op`<br>`double* y`<br>`double s1`<br>`double s2`<br>`uint32_t size` | `void` | Template helper that applies binary operation `op` with scale factor `s1` to even indices (`0, 2, 4, ...`) and `s2` to odd indices (`1, 3, 5, ...`) of array `y` up to `size`. |
| `dim` (Static) | `double* y`<br>`double s1`<br>`double s2`<br>`uint32_t size` | `void` | Converts an internal ODE array `y` to dimensional units by dividing even elements by `s1` and odd elements by `s2`. |
| `nonDim` (Static) | `double* y`<br>`double s1`<br>`double s2`<br>`uint32_t size` | `void` | Converts a dimensional ODE array `y` to non-dimensional units by multiplying even elements by `s1` and odd elements by `s2`. |

### `ParticleRestricted`
**Header:** `mantle_particle_restricted.hpp`

#### Purpose & Wrapped C Components
`ParticleRestricted` is an intermediate RAII wrapper constructed inside REBOUND callbacks (`m_additionalForces`, `m_collisionResolve`, `m_ode`, `m_heartbeat`) around a live C `reb_particle`'s Cartesian state and its `Meta` struct (`*static_cast<Meta*>(rp.ap)`). It manages converting the particle's state between REBOUND's internal non-dimensional units and physical dimensional units, and automatically re-nondimensionalizes any accessed ODE state arrays (`reb_ode.y`) in its destructor when the callback exits.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `state` | `m_cart_state_vector` | Copy of the particle's Cartesian state vector (extracted from `reb_particle` via `m_rebToCart`). |
| `meta_props` | `Meta&` | Reference to the particle's `Meta` instance stored in `reb_particle.ap`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `ParticleRestricted` (Constructor) | `m_cart_state_vector _state`<br>`Meta& m` | — | Initializes `state` with `_state`, binds `meta_props` to `m`, and sets the tracked ODE access counter `size` to `0`. |
| `operator=` | `const ParticleRestricted& other` | `ParticleRestricted&` | Copy-assignment operator that assigns `other.meta_props` to `meta_props` on non-self assignment. |
| `~ParticleRestricted` (Destructor) | None | — | Iterates over all ODE IDs recorded in `y_id` (accessed via `getOdey` during the lifetime of this instance) and calls `meta_props.nonDimId` to restore each `reb_ode.y` array to non-dimensional units. |
| `getParticleProperties` | None | `const std::any&` | Returns the type-erased custom particle properties from `meta_props`. |
| `getParticlePropertiesGen` | None | `const std::optional<PropertiesGen>&` | Returns the optional `PropertiesGen` struct from `meta_props`. |
| `getCollisionContext` | None | `const m_coll_context&` | Returns the collision context factory callback from `meta_props`. |
| `getAdditionalForces` | None | `const std::optional<std::vector<m_force_func>>&` | Returns the optional vector of additional force callbacks from `meta_props`. |
| `updateRebParticle` | None | `m_cart_state_vector` | Returns the current `state` vector so updated positions and velocities can be written back to `reb_particle`. |
| `getOdes` | None | `const std::optional<std::reference_wrapper<std::vector<Ode>>>` | Returns an optional reference wrapper to the particle's attached `Ode` vector. |
| `getOdey` | `uint32_t id` | `std::optional<double*>` | Retrieves and dimensionalizes the live `reb_ode.y` pointer for ODE `id` via `meta_props.getOdey(id)`, recording `id` in `y_id` so the destructor can restore non-dimensional scaling. |
| `getHash` | None | `uint32_t` | Returns the particle's hash (`meta_props.getHash()`). |
| `isRemovable` | None | `bool` | Returns whether the particle can be removed (`meta_props.isRemovable()`). |
| `incCollision` | None | `void` | Increments the particle's collision counter via `meta_props.incCollision()`. |
| `nonDim` | None | `void` | Converts `state` from dimensional to non-dimensional units by dividing by `meta_props.LU` and `meta_props.TU`. |
| `dim` | None | `void` | Converts `state` from non-dimensional to dimensional units by multiplying by `meta_props.LU` and `meta_props.TU`. |

### `ParticleLite`
**Header:** `mantle_particle_restricted.hpp`

#### Purpose & Wrapped C Components
`ParticleLite` is a lightweight facade over `ParticleRestricted` that is passed directly to user-defined callback functions (`m_force_func`, `m_coll_func`, `m_coll_context`, `m_ode_func`, and `m_remove_particle`). It exposes the particle's dimensionalized Cartesian state (`state`), type-safe accessors for custom properties, and coupled ODE state arrays while preventing user callbacks from mutating internal REBOUND pointers or simulation metadata directly.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `state` | `m_cart_state_vector&` | Reference to the underlying `ParticleRestricted::state` vector (in dimensional units during user callbacks). |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `ParticleLite` (Constructor) | `ParticleRestricted& p_` | — | Binds the internal reference `p` to `p_` and `state` to `p_.state`. |
| `setState` | `double x`, `double y`, `double z`<br>`double vx`, `double vy`, `double vz` | `void` | Updates the 3D position (`x, y, z`) and velocity (`vx, vy, vz`) components of `state` (used in collision resolution callbacks to modify post-collision kinematics). |
| `isType<T>` | None | `bool` | Template method that checks whether the stored `std::any` particle properties match type `T` (`typeid(T) == p.getParticleProperties().type()`). |
| `getParticleProperties<T>` | None | `std::optional<T>` | Template method that returns the custom particle properties cast to type `T` via `std::any_cast<T>` if `isType<T>()` is `true`, or `std::nullopt` otherwise. |
| `getParticlePropertiesGen` | None | `const std::optional<PropertiesGen>&` | Returns the optional `PropertiesGen` struct from the underlying `ParticleRestricted`. |
| `getOdey` | `uint32_t id` | `std::optional<double*>` | Returns a pointer to the dimensionalized live state array (`reb_ode.y`) for the coupled ODE matching `id`, or `std::nullopt` if not found. |

### `ParticleHistory`
**Header:** `mantle_particle_restricted.hpp`

#### Purpose & Wrapped C Components
`ParticleHistory` represents an immutable, dimensionalized snapshot of a single `reb_particle` and its `Meta` state at a specific output timestep during `Mantle::evalSim`. `Mantle::evalSim` returns a 2D vector (`m_history`, i.e., `std::vector<std::vector<ParticleHistory>>`) containing a vector of `ParticleHistory` objects for each recorded timestep.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `state` | `m_cart_state_vector` | The particle's Cartesian position, velocity, and simulation time (`x, y, z, vx, vy, vz, t`) at the snapshot timestep. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `ParticleHistory` (Default Constructor) | None | — | Default-constructs an empty `ParticleHistory` instance. |
| `ParticleHistory` (Constructor) | `m_cart_state_vector _state`<br>`Meta m` | — | Copies `_state` and metadata `m`, and immediately invokes `meta_props.syncOde()` to snapshot all coupled `reb_ode.y` arrays in dimensional units into `save_y`. |
| `collectHistory` | None | `bool` | Returns whether this particle is flagged for history collection (`meta_props.collectHistory()`). |
| `isType<T>` | None | `bool` | Template method that returns `true` if the stored custom particle properties match type `T`. |
| `getParticleProperties<T>` | None | `std::optional<T>` | Template method that returns the custom particle properties cast to `T` if `isType<T>()` is `true`, or `std::nullopt` otherwise. |
| `getParticlePropertiesGen` | None | `const std::optional<PropertiesGen>&` | Returns the optional `PropertiesGen` struct for this particle. |
| `getCollisionLog` | None | `int` | Returns the cumulative collision count for this particle up to this snapshot timestep. |
| `getOdey` | `uint32_t id` | `const std::optional<std::array<double, MAX_STATE>>` | Returns the saved dimensionalized ODE state vector (`save_y`) for the coupled ODE with identifier `id`, or `std::nullopt` if not found. |
| `getHash` | None | `uint32_t` | Returns the unique `reb_particle.hash` identifier of the particle. |
| `nonDim` | None | `void` | Scales `state` from dimensional to non-dimensional units using `meta_props.LU` and `meta_props.TU`. |
| `dim` | None | `void` | Scales `state` from non-dimensional to dimensional units using `meta_props.LU` and `meta_props.TU`. |

### `PropertiesGen`
**Header:** `mantle_particle_types.hpp`

#### Purpose & Wrapped C Components
`PropertiesGen` is a plain data structure providing a standardized set of physical, geometric, thermal, and optical particle parameters. It is stored inside `Meta` (attached to `reb_particle.ap`) and `Ode` (attached to `reb_ode.ref`) so custom radiation, drag, or thermal ODE callbacks can access common physical properties without requiring `std::any_cast`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `mass` | `double` | Particle mass. |
| `radius` | `double` | Particle physical radius. |
| `Ac` | `double` | Cross-sectional area. |
| `As` | `double` | Surface area. |
| `B` | `double` | Ballistic or radiation pressure coefficient parameter. |
| `Cp` | `double` | Specific heat capacity. |
| `cold_transparency` | `double` | Optical transparency in the cold state (below `transition_temp`). |
| `cold_absorptivity` | `double` | Optical absorptivity in the cold state. |
| `cold_reflectivity` | `double` | Optical reflectivity in the cold state. |
| `cold_emissivity` | `double` | Thermal emissivity in the cold state. |
| `hot_transparency` | `double` | Optical transparency in the hot state (above `transition_temp`). |
| `hot_absorptivity` | `double` | Optical absorptivity in the hot state. |
| `hot_reflectivity` | `double` | Optical reflectivity in the hot state. |
| `hot_emissivity` | `double` | Thermal emissivity in the hot state. |
| `transition_temp` | `double` | Phase or state transition temperature threshold separating cold and hot optical regimes. |

#### Public Methods
None (aggregate struct).

### `m_cart_state_vector`
**Header:** `mantle_particle.hpp`

#### Purpose & Wrapped C Components
Represents a 3D Cartesian position and velocity state vector at time `t`. It maps directly to the Cartesian kinematic fields (`x, y, z, vx, vy, vz`) of REBOUND's `reb_particle` struct combined with the simulation time `r->t`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `x` | `double` | Position component along the X-axis. |
| `y` | `double` | Position component along the Y-axis. |
| `z` | `double` | Position component along the Z-axis. |
| `vx` | `double` | Velocity component along the X-axis. |
| `vy` | `double` | Velocity component along the Y-axis. |
| `vz` | `double` | Velocity component along the Z-axis. |
| `t` | `double` | Simulation time associated with the state. |

#### Public Methods
None (aggregate struct).

### `m_kep_state_vector`
**Header:** `mantle_particle.hpp`

#### Purpose & Wrapped C Components
Represents a set of classical Keplerian orbital elements at time `t`. When a `Particle` is initialized with `m_kep_state_vector`, `Mantle` passes these elements to REBOUND's `reb_particle_from_orbit(r->G, r->particles[0], m, a, e, inc, omega, Omega, f)` function to compute the corresponding `reb_particle` Cartesian state around the primary body (`r->particles[0]`).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `a` | `double` | Semi-major axis. |
| `e` | `double` | Orbital eccentricity. |
| `inc` | `double` | Orbital inclination (radians). |
| `Omega` | `double` | Longitude of the ascending node $\Omega$ (radians). |
| `omega` | `double` | Argument of periapsis $\omega$ (radians). |
| `f` | `double` | True anomaly $f$ (radians). |
| `t` | `double` | Epoch or simulation time associated with the orbital elements. |

#### Public Methods
None (aggregate struct).

### `m_accel`
**Header:** `mantle_particle.hpp`

#### Purpose & Wrapped C Components
Represents a 3D acceleration vector returned by user-defined force callbacks (`m_force_func`). Within `m_additionalForces` and `m_additionalForcesPar`, `m_accel` vectors from all force functions on a particle are summed, non-dimensionalized, and added to `reb_particle.ax`, `reb_particle.ay`, and `reb_particle.az`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `ax` | `double` | Acceleration component along the X-axis. |
| `ay` | `double` | Acceleration component along the Y-axis. |
| `az` | `double` | Acceleration component along the Z-axis. |

#### Public Methods & Friend Operators

| Operator | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `operator+` (Friend) | `const m_accel& t`<br>`const m_accel& u` | `m_accel` | Returns the component-wise vector sum `{t.ax + u.ax, t.ay + u.ay, t.az + u.az}`. |
| `operator*` (Friend) | `double s`<br>`const m_accel& u` | `m_accel` | Scales the acceleration vector by scalar `s` (`{s * u.ax, s * u.ay, s * u.az}`). |
| `operator*` (Friend) | `const m_accel& u`<br>`double s` | `m_accel` | Commutative scalar multiplication returning `{u.ax * s, u.ay * s, u.az * s}`. |

### `Scale<Op>`
**Header:** `mantle_particle.hpp`

#### Purpose & Wrapped C Components
A templated functor class (compatible with `std::visit` on `Particle::state`) that converts `m_cart_state_vector` and `m_kep_state_vector` objects between physical dimensional units and REBOUND's internal non-dimensional units using binary operator `Op` (`std::divides<double>` for non-dimensionalization, `std::multiplies<double>` for dimensionalization).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `op` | `Op` | Binary arithmetic functor (`std::divides<double>` or `std::multiplies<double>`). |
| `LU` | `double` | Characteristic length unit. |
| `TU` | `double` | Characteristic time unit. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `Scale` (Constructor) | `Op o`<br>`double l`<br>`double t` | — | Initializes the scaling functor with binary operator `o`, length unit `LU = l`, and time unit `TU = t`. |
| `operator()` | `m_cart_state_vector& state` | `void` | Scales Cartesian positions (`x, y, z`) by `LU`, velocities (`vx, vy, vz`) by `LU / TU`, and time (`t`) by `TU` in-place using `op`. |
| `operator()` | `m_kep_state_vector& state` | `void` | Scales Keplerian semi-major axis (`a`) by `LU` and time (`t`) by `TU` in-place using `op` (angular elements and eccentricity remain dimensionless). |

### `VisitState`
**Source:** `mantle.cpp`

#### Purpose & Wrapped C Components
Internal `std::visit` visitor struct used by `Mantle::addParticles()` to initialize and register a C `reb_particle` in `reb_simulation` from a `Particle`'s `std::variant<m_cart_state_vector, m_kep_state_vector>` state. Wraps `reb_simulation_add` and `reb_particle_from_orbit`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `r` | `reb_simulation*` | Pointer to the target REBOUND simulation. |
| `m` | `double` | Non-dimensionalized mass of the particle. |
| `rp` | `double` | Non-dimensionalized physical radius of the particle. |
| `meta` | `const Meta&` | Reference to the particle's `Meta` object, whose address is stored in `reb_particle.ap`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `operator()` | `m_cart_state_vector state` | `void` | Populates a `reb_particle` directly with Cartesian coordinates (`x, y, z, vx, vy, vz`), mass `m`, radius `rp`, `meta.hash`, and `ap = &meta`, then adds it via `reb_simulation_add(r, p)`. |
| `operator()` | `m_kep_state_vector state` | `void` | Generates a `reb_particle` around `r->particles[0]` via `reb_particle_from_orbit(r->G, r->particles[0], m, state.a, state.e, state.inc, state.omega, state.Omega, state.f)`, sets `r = rp`, `hash = meta.hash`, and `ap = &meta`, then adds it via `reb_simulation_add(r, p)`. |

## Simulation Configuration Classes

All configuration classes inherit from `MConfig` and are passed as `std::vector<std::unique_ptr<MConfig>>` to the `Mantle` constructor to configure the `Mantle` wrapper and the underlying C `reb_simulation` struct.

### `MConfig`
**Header:** `mantle_config.hpp`

#### Purpose & Wrapped C Components
Abstract base polymorphic struct for all simulation configuration modules. Each derived struct implements `updateSim` to mutate `Mantle` and the C `reb_simulation` struct.

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` (Pure Virtual) | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Applies configuration settings to the `Mantle` wrapper `m` and the REBOUND `reb_simulation` pointer `r`. Returns an empty string `""` on success or an error message string. |
| `~MConfig` (Virtual Destructor) | None | — | Default virtual destructor for safe polymorphic cleanup. |

### `General`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Configures global simulation parameters, characteristic scaling units, callback hooks, and OpenMP multithreading. It wraps general fields of `reb_simulation` (`r->G`, `r->testparticle_hidewarnings`, `r->dt`, `r->N_active`, `r->softening`) as well as OpenMP's `omp_set_num_threads` when compiled with `OPENMP`.

#### Public Data Members

| Member | Type | Default Value | Description |
| :--- | :--- | :--- | :--- |
| `G` | `std::optional<double>` | `std::nullopt` | Gravitational constant (`r->G`). |
| `MU` | `double` | `1.0` | Characteristic mass unit stored on `Mantle::MU`. |
| `LU` | `double` | `1.0` | Characteristic length unit stored on `Mantle::LU`. |
| `TU` | `double` | `1.0` | Characteristic time unit stored on `Mantle::TU`. |
| `silence` | `bool` | `true` | Suppresses REBOUND test-particle warning messages (`r->testparticle_hidewarnings`). |
| `has_ode` | `bool` | `false` | Enables custom coupled ODE integration (`Mantle::has_ode`). |
| `has_additional_forces` | `bool` | `false` | Enables custom non-gravitational force evaluation (`Mantle::has_additional_forces`). |
| `force_is_velocity_dep` | `bool` | `true` | Sets whether additional forces depend on particle velocity (`r->force_is_velocity_dependent`). |
| `dt` | `std::optional<double>` | `std::nullopt` | Initial integration timestep (`r->dt`). |
| `active` | `std::optional<int>` | `std::nullopt` | Number of active massive particles (`r->N_active`); remaining particles act as test particles. |
| `softening` | `std::optional<double>` | `std::nullopt` | Gravitational softening length (`r->softening`). |
| `threads` | `std::optional<int>` | `std::nullopt` | Number of OpenMP threads to configure via `omp_set_num_threads` (when compiled with `OPENMP`). |
| `removeParticleTest` | `std::optional<m_remove_particle>` | `std::nullopt` | Heartbeat callback predicate used to test and remove particles during integration. |
| `simTerminate` | `std::optional<m_sim_terminate>` | `std::nullopt` | Heartbeat callback predicate used to check whether the simulation should terminate early. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `General` (Constructor) | None | — | Initializes default flags (`has_ode = false`, `has_additional_forces = false`, `force_is_velocity_dep = true`, `silence = true`) and unit scaling factors (`MU = 1.0`, `LU = 1.0`, `TU = 1.0`). |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Applies all populated general settings, unit scales, and heartbeat callbacks to `m` and `r`, and sets the OpenMP thread count if `OPENMP` is defined. Returns `""`. |

### `SimBox`
**Header:** `mantle_config.hpp`

#### Purpose & Wrapped C Components
Helper data structure that stores spatial domain parameters for octree-based gravity and collision algorithms (`TreeGravity`, `TreeCollision`, `LineTreeCollision`). Its fields are passed directly to REBOUND's `reb_simulation_configure_box(r, box_size, x_boxes, y_boxes, z_boxes)`.

#### Public Data Members

| Member | Type | Default Value | Description |
| :--- | :--- | :--- | :--- |
| `box_size` | `std::optional<double>` | `std::nullopt` | Side length of a single root simulation box. |
| `x_boxes` | `std::optional<int>` | `0` | Number of root boxes along the X-axis. |
| `y_boxes` | `std::optional<int>` | `0` | Number of root boxes along the Y-axis. |
| `z_boxes` | `std::optional<int>` | `0` | Number of root boxes along the Z-axis. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `SimBox` (Constructor) | None | — | Initializes `x_boxes = 0`, `y_boxes = 0`, and `z_boxes = 0`, leaving `box_size` unset (`std::nullopt`). |

## Gravity Configuration Classes

### `BasicGravity`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures the gravity solver mode on `reb_simulation.gravity`. Note that in `mantle_config.cpp`, `BasicGravity::updateSim` sets `r->gravity = reb_simulation::REB_GRAVITY_TREE`.

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->gravity = reb_simulation::REB_GRAVITY_TREE` and returns `""`. |

### `CompensatedGravity`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's Kahan-compensated direct summation gravity solver (`reb_simulation::REB_GRAVITY_COMPENSATED`).

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->gravity = reb_simulation::REB_GRAVITY_COMPENSATED` and returns `""`. |

### `TreeGravity`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's Barnes-Hut octree gravity solver (`reb_simulation::REB_GRAVITY_TREE`). It wraps `r->opening_angle2`, calls `reb_simulation_configure_box`, and calls `reb_mpi_init(r)` when compiled with `MPI`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `box` | `SimBox` | Simulation domain box configuration passed to `reb_simulation_configure_box`. |
| `opening_angle2` | `std::optional<double>` | Square of the cell opening angle $\theta^2$ (`r->opening_angle2`). |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->gravity = reb_simulation::REB_GRAVITY_TREE`, sets `r->opening_angle2` if provided, validates `box.box_size` (returning `"Error: No box size set"` if missing), calls `reb_simulation_configure_box`, initializes MPI if `MPI` is defined, and returns `""`. |

### `JacobiGravity`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's Jacobi coordinate gravity solver (`reb_simulation::REB_GRAVITY_JACOBI`), typically used alongside symplectic integrators.

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->gravity = reb_simulation::REB_GRAVITY_JACOBI` and returns `""`. |

### `NoGravity`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and disables gravitational interactions in REBOUND (`reb_simulation::REB_GRAVITY_NONE`).

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->gravity = reb_simulation::REB_GRAVITY_NONE` and returns `""`. |

## Collision Detection Configuration Classes

### `DirectCollision`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's direct $O(N^2)$ collision detection (`reb_simulation::REB_COLLISION_DIRECT`). It sets the default collision resolver to `reb_collision_resolve_hardsphere`, or delegates to `m_collisionResolve` if a custom `collisionContext` is provided.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `collisionContext` | `std::optional<m_coll_context>` | Optional custom collision resolution factory function stored on `Mantle::collisionContext`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->collision = reb_simulation::REB_COLLISION_DIRECT`, sets `r->collision_resolve = reb_collision_resolve_hardsphere`, copies `collisionContext` to `m->collisionContext` if present, and returns `""`. |

### `LineCollision`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's continuous line-trajectory collision detection (`reb_simulation::REB_COLLISION_LINE`) along with `reb_collision_resolve_hardsphere` (or a custom `collisionContext`).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `collisionContext` | `std::optional<m_coll_context>` | Optional custom collision resolution factory function stored on `Mantle::collisionContext`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->collision = reb_simulation::REB_COLLISION_LINE`, sets `r->collision_resolve = reb_collision_resolve_hardsphere`, copies `collisionContext` to `m->collisionContext` if present, and returns `""`. |

### `TreeCollision`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's octree collision search (`reb_simulation::REB_COLLISION_TREE`). Wraps `reb_collision_resolve_hardsphere`, `reb_simulation_configure_box`, and `reb_mpi_init` (when `MPI` is enabled).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `box` | `SimBox` | Simulation box dimensions passed to `reb_simulation_configure_box`. |
| `collisionContext` | `std::optional<m_coll_context>` | Optional custom collision resolution factory function stored on `Mantle::collisionContext`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->collision = reb_simulation::REB_COLLISION_TREE` and `r->collision_resolve = reb_collision_resolve_hardsphere`, assigns `m->collisionContext` if set, validates `box.box_size` (returning `"Error: No box size set"` if absent), calls `reb_simulation_configure_box`, calls `reb_mpi_init(r)` if `MPI` is enabled, and returns `""`. |

### `LineTreeCollision`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's octree continuous line-trajectory collision detection (`reb_simulation::REB_COLLISION_LINETREE`). Wraps `reb_collision_resolve_hardsphere`, `reb_simulation_configure_box`, and `reb_mpi_init` (when `MPI` is enabled).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `box` | `SimBox` | Simulation box dimensions passed to `reb_simulation_configure_box`. |
| `collisionContext` | `std::optional<m_coll_context>` | Optional custom collision resolution factory function stored on `Mantle::collisionContext`. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->collision = reb_simulation::REB_COLLISION_LINETREE` and `r->collision_resolve = reb_collision_resolve_hardsphere`, assigns `m->collisionContext` if set, validates `box.box_size` (returning `"Error: No box size set"` if absent), calls `reb_simulation_configure_box`, calls `reb_mpi_init(r)` if `MPI` is enabled, and returns `""`. |

### `NoCollision`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and disables collision detection in REBOUND (`reb_simulation::REB_COLLISION_NONE`).

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->collision = reb_simulation::REB_COLLISION_NONE` and returns `""`. |

## Integrator Configuration Classes

### `IAS15Integrator`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's 15th-order adaptive Gauss-Radau integrator (`reb_simulation::REB_INTEGRATOR_IAS15`) and wraps the `reb_integrator_ias15` sub-struct (`r->ri_ias15`).

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `epsilon` | `std::optional<double>` | Relative error tolerance parameter (`r->ri_ias15.epsilon`). |
| `min_dt` | `std::optional<double>` | Minimum allowed adaptive timestep (`r->ri_ias15.min_dt`). |
| `adaptive_mode` | `std::optional<bool>` | Adaptive timestepping mode selector (`r->ri_ias15.adaptive_mode`). |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->integrator = reb_simulation::REB_INTEGRATOR_IAS15` and copies `epsilon`, `min_dt`, and `adaptive_mode` into `r->ri_ias15` if present. Returns `""`. |

### `BSIntegrator`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's Gragg-Bulirsch-Stoer integrator (`reb_simulation::REB_INTEGRATOR_BS`) and wraps the `reb_integrator_bs` sub-struct (`r->ri_bs`). Can be configured either as the full N-body + ODE integrator (`ode_only = false`) or exclusively as the sub-integrator for coupled `reb_ode` systems (`ode_only = true`, the default).

#### Public Data Members

| Member | Type | Default Value | Description |
| :--- | :--- | :--- | :--- |
| `ode_only` | `bool` | `true` | If `false`, sets `r->integrator = reb_simulation::REB_INTEGRATOR_BS` for the entire simulation; if `true`, configures `r->ri_bs` tolerances for coupled ODEs while leaving the primary N-body integrator unchanged. |
| `eps_rel` | `std::optional<double>` | `std::nullopt` | Relative error tolerance (`r->ri_bs.eps_rel`). |
| `eps_abs` | `std::optional<double>` | `std::nullopt` | Absolute error tolerance (`r->ri_bs.eps_abs`). |
| `min_dt` | `std::optional<double>` | `std::nullopt` | Minimum allowed step size (`r->ri_bs.min_dt`). |
| `max_dt` | `std::optional<double>` | `std::nullopt` | Maximum allowed step size (`r->ri_bs.max_dt`). |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `BSIntegrator` (Constructor) | None | — | Initializes `ode_only = true`. |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->integrator = reb_simulation::REB_INTEGRATOR_BS` when `ode_only` is `false`, and updates `r->ri_bs.eps_rel`, `eps_abs`, `max_dt`, and `min_dt` if provided. Returns `""`. |

### `WHFastIntegrator`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and selects REBOUND's Wisdom-Holman symplectic integrator (`reb_simulation::REB_INTEGRATOR_WHFAST`). Note: While the struct declares optional fields corresponding to `r->ri_whfast` settings, `WHFastIntegrator::updateSim` in `mantle_config.cpp` currently sets `r->integrator = reb_simulation::REB_INTEGRATOR_WHFAST` and returns `""`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `corrector` | `std::optional<uint32_t>` | Order of the first symplectic corrector. |
| `corrector2` | `std::optional<bool>` | Enables or disables the second symplectic corrector. |
| `kernel` | `std::optional<WHFastKernel>` | WHFast kernel method (`WHFastKernel` enum). |
| `coordinates` | `std::optional<WHFastCoord>` | Coordinate system used by WHFast (`WHFastCoord` enum). |
| `recalculate_coordinates_this_timestep` | `std::optional<bool>` | Flag to trigger coordinate recalculation on the current step. |
| `safe_mode` | `std::optional<bool>` | Enables or disables WHFast safe mode (synchronizing every timestep). |
| `keep_unsynchronized` | `std::optional<bool>` | Flag to keep particles unsynchronized after outputs. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->integrator = reb_simulation::REB_INTEGRATOR_WHFAST` and returns `""`. |

### `LeapFrogIntegrator`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and configures REBOUND's second-order symplectic Leapfrog integrator (`reb_simulation::REB_INTEGRATOR_LEAPFROG`).

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG` and returns `""`. |

### `SEIIntegrator`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and selects REBOUND's Symplectic Epicycle Integrator for shearing-sheet simulations (`reb_simulation::REB_INTEGRATOR_SEI`). Note: While `omega` and `omegaz` are declared as public members, `SEIIntegrator::updateSim` in `mantle_config.cpp` currently sets `r->integrator = reb_simulation::REB_INTEGRATOR_SEI` and returns `""`.

#### Public Data Members

| Member | Type | Description |
| :--- | :--- | :--- |
| `omega` | `std::optional<double>` | Epicyclic or orbital frequency $\Omega$. |
| `omegaz` | `std::optional<double>` | Vertical epicyclic frequency $\Omega_z$. |

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->integrator = reb_simulation::REB_INTEGRATOR_SEI` and returns `""`. |

### `NoIntegrator`
**Header:** `mantle_config.hpp` | **Source:** `mantle_config.cpp`

#### Purpose & Wrapped C Components
Inherits from `MConfig` and disables particle motion integration (`reb_simulation::REB_INTEGRATOR_NONE`).

#### Public Data Members
None.

#### Public Methods

| Method | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `updateSim` | `class Mantle* m`<br>`reb_simulation* r` | `std::string` | Sets `r->integrator = reb_simulation::REB_INTEGRATOR_NONE` and returns `""`. |

## Supporting Enumerations, Type Aliases & Free C-Bridge Functions

### Enumerations & Constants

| Identifier | Kind | Values / Definition | Description |
| :--- | :--- | :--- | :--- |
| `WHFastKernel` | `enum class` | `REB_WHFAST_KERNEL_DEFAULT`, `REB_WHFAST_KERNEL_MODIFIEDKICK`, `REB_WHFAST_KERNEL_COMPOSITION`, `REB_WHFAST_KERNEL_LAZY` | Kernel options for `WHFastIntegrator` (`mantle_config.hpp`). |
| `WHFastCoord` | `enum class` | `REB_WHFAST_COORDINATES_JACOBI`, `REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC`, `REB_WHFAST_COORDINATES_WHDS`, `REB_WHFAST_COORDINATES_BARYCENTRIC` | Coordinate system options for `WHFastIntegrator` (`mantle_config.hpp`). |
| `MAX_STATE` | Macro | `10` | Maximum number of state variables supported per `Ode` instance (`mantle_particle.hpp`). |
| `MAX_ODE_CALL` | Macro | `10` | Maximum number of `getOdey` calls tracked per `ParticleRestricted` instance (`mantle_particle.hpp`). |

### Type Aliases & Callback Signatures

| Alias | Underlying Type | Description |
| :--- | :--- | :--- |
| `m_history` | `std::vector<std::vector<ParticleHistory>>` | 2D vector of `ParticleHistory` snapshots indexed by `[timestep][particle]` returned by `Mantle::evalSim`. |
| `m_force_func` | `m_accel (*)(const ParticleLite&)` | Callback signature for computing an additional dimensional acceleration `m_accel` on a particle. |
| `m_coll_func` | `int (*)(ParticleLite&, ParticleLite&)` | Callback signature for resolving a collision between two particles; returns the REBOUND collision resolution outcome code (`0`, `1`, `2`, or `3`). |
| `m_coll_context` | `m_coll_func (*)(const ParticleLite&, const ParticleLite&)` | Higher-order dispatcher callback that inspects two colliding `ParticleLite` objects and returns the appropriate `m_coll_func` resolver. |
| `m_ode_func` | `void (*)(double* const yDot, const double* const y, double t, const ParticleLite&)` | Callback signature for computing dimensional ODE state derivatives `yDot` given current dimensional state `y`, dimensional time `t`, and a `ParticleLite` reference. |
| `m_remove_particle` | `bool (*)(const ParticleLite&)` | Predicate callback evaluated at each heartbeat step; returns `true` if the particle should be removed from the simulation. |
| `m_sim_terminate` | `bool (*)(int)` | Predicate callback evaluated at each heartbeat step given the active particle count `r->N`; returns `true` to stop the simulation early. |

### Internal C-Callback Bridge Functions (`mantle.hpp` / `mantle.cpp`)

| Function | Parameters | Return Type | Description |
| :--- | :--- | :--- | :--- |
| `m_additionalForces` | `struct reb_simulation* const r` | `void` | Sequential C callback bound to `r->additional_forces` (when `OPENMP` is not defined). Dimensionalizes each particle, accumulates accelerations across all attached `m_force_func` callbacks, converts the total acceleration back to non-dimensional units (`TU * TU / LU`), and adds it to `rp.ax`, `rp.ay`, `rp.az`. |
| `m_additionalForcesPar` | `struct reb_simulation* const r` | `void` | OpenMP-parallelized version (`#pragma omp parallel for`) of `m_additionalForces`, bound when `OPENMP` is defined. |
| `m_collisionResolve` | `struct reb_simulation* r`<br>`struct reb_collision c` | `int` | C callback bound to `r->collision_resolve`. Wraps colliding particles `c.p1` and `c.p2` in `ParticleRestricted` and `ParticleLite`, invokes `collisionContext` and the returned `m_coll_func`, writes updated states back to `reb_particle`, increments collision counters (`p1.incCollision()`, `p2.incCollision()`, `r->collisions_log_n++`), and returns the REBOUND removal code. |
| `m_ode` | `struct reb_ode* const ode`<br>`double* const yDot`<br>`const double* const y`<br>`const double t` | `void` | C callback bound to `ho->derivatives`. Looks up the host `reb_particle` via `reb_simulation_particle_by_hash`, dimensionalizes `y` and the particle state, invokes the user `m_ode_func`, and non-dimensionalizes both `yDot` and `y`. |
| `m_heartbeat` | `struct reb_simulation* r` | `void` | Sequential C callback bound to `r->heartbeat`. Evaluates `removeParticleTest` on each removable particle, frees any attached `reb_ode` instances and calls `reb_simulation_remove_particle_by_hash` on matches, and calls `reb_simulation_stop(r)` if `simTerminate(r->N)` returns `true`. |
| `m_heartbeatPar` | `struct reb_simulation* r` | `void` | OpenMP-parallelized version of `m_heartbeat` that gathers removable particles across threads before removing them and checking `simTerminate(r->N)`. |
| `m_rebToCart` | `struct reb_particle& p`, `double t_`<br>*(or `struct reb_particle* p`, `double t_`)* | `m_cart_state_vector` | Utility functions that extract `{x, y, z, vx, vy, vz, t_}` from a C `reb_particle` reference or pointer into an `m_cart_state_vector`. |
