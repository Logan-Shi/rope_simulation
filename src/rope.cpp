#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass,
           float k, vector<int> pinned_nodes) {
  // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and
  // containing `num_nodes` nodes.
  //
  // New masses can be created (see mass.h):
  // Mass *m = new Mass(position, mass, bool)
  // and should be added to the masses vector.
  //
  // Springs can be created (see spring.h):
  // Spring *s = new Spring(mass_a, mass_b, spring_const)
  // and should be added to the springs vector.
  //
  // Masses corresponding to indices in pinned_nodes
  // should have their pinned field set to true.

  Vector2D incre = (start - end)/(num_nodes - 1);
  for (int i = 0; i < num_nodes; ++i)
  {
    Vector2D node_pos = start + i*incre;
    bool is_pinned = std::find(pinned_nodes.begin(), pinned_nodes.end(), i) != pinned_nodes.end();
    Mass *m = new Mass(node_pos, node_mass, is_pinned);
    masses.push_back(m);
  }

  for (int i = 0; i < num_nodes - 1; ++i)
  {
    Spring *s = new Spring(masses[i], masses[i + 1], k);
    springs.push_back(s);
  }
}

void Rope::simulateEuler(float delta_t, Vector2D gravity) {
  for (auto &s : springs) {
    // std::cout<<"m1 position x: "<<s->m1->position[0]<<"\n";
    // std::cout<<"m2 position x: "<<s->m2->position[0]<<"\n";
    // TODO (Part 2.1): Use Hooke's law to calculate the force on a node

    Vector2D m1tom2 = s->m2->position - s->m1->position;
    double length = m1tom2.norm();
    Vector2D node_forces = -s->k * m1tom2 / length * (length - s->rest_length);

    // TODO (Part 4.1): Add damping forces
    s->m1->forces += -node_forces;
    s->m2->forces += node_forces;

    Vector2D relative_vel = s->m2->velocity - s->m1->velocity;
    Vector2D damping = 0.5*relative_vel;
    s->m1->forces += damping;
    s->m2->forces += -damping;
  }

  for (auto &m : masses) {
    // std::cout<<"position x: "<<m->position[0]<<"\n";
    // std::cout<<"forces x: "<<m->forces[0]<<"\n";
    if (!m->pinned) {
      // TODO (Part 2.1): Add the force due to gravity, then compute the new
      // velocity and position

      // m->velocity = m->last_position - m->position;
      // m->last_position = m->position;

      m->velocity = m->velocity + (m->forces/m->mass + gravity) * delta_t;
      m->position = m->position + m->velocity * delta_t;
    }

    // TODO Reset all forces on each mass
    m->forces = Vector2D(0,0);
  }
}

void Rope::simulateVerlet(float delta_t, Vector2D gravity) {
  // TODO (Part 3.1): Clear forces
  for (auto &s : springs) {
    // TODO (Part 3.1): Simulate one timestep of the rope using explicit Verlet
    // std::cout<<"m1 position x: "<<s->m1->position[0]<<"\n";
    // std::cout<<"m2 position x: "<<s->m2->position[0]<<"\n";
    // TODO (Part 2.1): Use Hooke's law to calculate the force on a node

    Vector2D m1tom2 = s->m2->position - s->m1->position;
    double length = m1tom2.norm();
    Vector2D node_forces = -s->k * m1tom2 / length * (length - s->rest_length);

    // TODO (Part 4.1): Add damping forces
    s->m1->forces += -node_forces;
    s->m2->forces += node_forces;
  }

  for (auto &m : masses) {
    if (!m->pinned) {
      Vector2D temp_position = m->position;

      // TODO (Part 3.1): Set the new position of the rope mass
      // TODO (Part 4.2): Add global Verlet damping

      m->position = m->position + 0.9995 * (m->position - m->last_position) + (m->forces/m->mass + gravity) * delta_t * delta_t;
      m->last_position = temp_position;
    }
    m->forces = Vector2D(0,0);
  }
}
}
