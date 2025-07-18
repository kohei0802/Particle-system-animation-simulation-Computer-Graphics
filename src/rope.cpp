#include <iostream>
#include <vector>

#include "CGL/vector2D.h"
#include <assert.h>

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.

        assert(num_nodes > 1);

        masses = vector<Mass*>();
        masses.resize(num_nodes);
        springs = vector<Spring*>();
        springs.resize(num_nodes - 1);

        for (int i = 0; i < num_nodes; i++) 
        {
            Vector2D position = start + (end - start)  * ((double) i / (double) (num_nodes - 1));

            masses[i] = new Mass(position, node_mass, false);
            masses[i]->forces = Vector2D(0, 0);
            masses[i]->velocity = Vector2D(0, 0);
            masses[i]->last_position = position;
        }

       for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
       }

       for (int i = 0; i < num_nodes - 1; i++) 
       {
            springs[i] = new Spring(masses[i], masses[i+1], k);
       }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        // 1.2

        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            double k = s->k;
            double l = s->rest_length;
            Vector2D direction = s->m2->position - s->m1->position;
            double norm = direction.norm();
        
            Vector2D dir_normed = direction / norm;
            Vector2D relative_vel = s->m2->velocity - s->m1->velocity;


            Vector2D f12 = k * direction * (norm - l) / norm;
            Vector2D f21 = -k * direction * (norm - l) / norm;
        

            s->m1->forces += f12;
            s->m2->forces += f21;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position

                // implicit Euler
                m->forces += gravity;

                Vector2D accel = m->forces / m->mass;
                m->velocity = m->velocity + accel * delta_t;
                m->position = m->position + m->velocity * delta_t;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        // 1.3, 1.4

        double damping_factor = 5e-5;

        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            double k = s->k;
            double l = s->rest_length;
            Vector2D direction = s->m2->position - s->m1->position;
            double norm = direction.norm();
        
            Vector2D dir_normed = direction / norm;
            Vector2D relative_vel = s->m2->velocity - s->m1->velocity;


            Vector2D f12 = k * direction * (norm - l) / norm;
            Vector2D f21 = -k * direction * (norm - l) / norm;
        

            s->m1->forces += f12;
            s->m2->forces += f21;
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass

                m->forces += gravity;

                Vector2D accel = m->forces / m->mass;

                m->position = temp_position + (1 - damping_factor) * (temp_position - m->last_position) + accel * delta_t * delta_t;

                m->last_position = temp_position;
                
                // TODO (Part 4): Add global Verlet damping

                m->forces = Vector2D(0, 0);
            }
        }
    }
}
