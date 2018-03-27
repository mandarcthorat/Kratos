//   $Main author: Guillermo Casas
//

// Project includes

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_particle_watcher.h"

namespace Kratos
{

typedef ModelPart::ElementsContainerType::iterator ElementsIteratorType;
typedef Kratos::AnalyticSphericParticle AnalyticParticle;

void AnalyticParticleWatcher::MakeMeasurements(ModelPart& analytic_model_part)
{
    const double current_time = analytic_model_part.GetProcessInfo()[TIME];    

    std::vector<InterParticleImpactDataOfAllParticlesSingleTimeStep> time_step_database_by_threads(OpenMPUtils::GetNumThreads(), current_time);
    std::vector<FaceParticleImpactDataOfAllParticlesSingleTimeStep> face_time_step_database_by_threads(OpenMPUtils::GetNumThreads(), current_time);

    std::vector<std::map<int, InterParticleImpactDataOfAllTimeStepsSingleParticle> > inter_particle_impact_data_of_all_time_steps_by_threads;
    inter_particle_impact_data_of_all_time_steps_by_threads.resize(OpenMPUtils::GetNumThreads());

    std::vector<std::map<int, FaceParticleImpactDataOfAllTimeStepsSingleParticle> > face_particle_impact_data_of_all_time_steps_by_threads;
    face_particle_impact_data_of_all_time_steps_by_threads.resize(OpenMPUtils::GetNumThreads());    

    #pragma omp parallel for if(analytic_model_part.NumberOfElements(0)>100)
    for (int k=0; k<(int)analytic_model_part.NumberOfElements(0); k++) {
        ElementsIteratorType i_elem = analytic_model_part.ElementsBegin() + k;
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        const int n_collisions = particle.GetNumberOfCollisions();

        const int id = int(i_elem->Id());
        InterParticleImpactDataOfAllTimeStepsSingleParticle& particle_database = GetParticleDataBase(id, inter_particle_impact_data_of_all_time_steps_by_threads[OpenMPUtils::ThisThread()]);
        if (n_collisions){            
            const array_1d<int, 4> &colliding_ids = particle.GetCollidingIds();
            const array_1d<double, 4> &colliding_normal_vel = particle.GetCollidingNormalRelativeVelocity();
            const array_1d<double, 4> &colliding_tangential_vel = particle.GetCollidingTangentialRelativeVelocity();
            const array_1d<double, 4> &colliding_linear_impulse = particle.GetCollidingLinearImpulse();

            for (int i = 0; i < n_collisions; ++i){
                time_step_database_by_threads[OpenMPUtils::ThisThread()].PushBackImpacts(id, colliding_ids[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
                particle_database.PushBackImpacts(current_time, colliding_ids[i], colliding_normal_vel[i], colliding_tangential_vel[i], colliding_linear_impulse[i]);
            }
        }

        const int n_collisions_with_walls = particle.GetNumberOfCollisionsWithFaces();
        FaceParticleImpactDataOfAllTimeStepsSingleParticle& flat_wall_particle_database = GetParticleFaceDataBase(id, face_particle_impact_data_of_all_time_steps_by_threads[OpenMPUtils::ThisThread()]);
        if (n_collisions_with_walls){                       
            const array_1d<int, 4> &colliding_ids_with_walls = particle.GetCollidingFaceIds();
            const array_1d<double, 4> &colliding_normal_vel = particle.GetCollidingFaceNormalRelativeVelocity();
            const array_1d<double, 4> &colliding_tangential_vel = particle.GetCollidingFaceTangentialRelativeVelocity();

            for (int i = 0; i < n_collisions_with_walls; ++i){
                face_time_step_database_by_threads[OpenMPUtils::ThisThread()].PushBackImpacts(id, colliding_ids_with_walls[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
                flat_wall_particle_database.PushBackImpacts(current_time, colliding_ids_with_walls[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
            }
        }
    }


    //Now we collect the info from the threads:
    InterParticleImpactDataOfAllParticlesSingleTimeStep time_step_database(current_time);
    FaceParticleImpactDataOfAllParticlesSingleTimeStep face_time_step_database(current_time);

    for (int i=0; i<OpenMPUtils::GetNumThreads(); i++) {
        time_step_database.PushBackImpacts(time_step_database_by_threads[i]);
        face_time_step_database.PushBackImpacts(face_time_step_database_by_threads[i]);
        for(std::map<int, InterParticleImpactDataOfAllTimeStepsSingleParticle>::iterator it = inter_particle_impact_data_of_all_time_steps_by_threads[i].begin(); it != inter_particle_impact_data_of_all_time_steps_by_threads[i].end(); ++it) {
            if (mInterParticleImpactDataOfAllTimeSteps.find(it->first) == mInterParticleImpactDataOfAllTimeSteps.end() ) {
                mInterParticleImpactDataOfAllTimeSteps[it->first] = it->second;
            }             
        }
        for(std::map<int, FaceParticleImpactDataOfAllTimeStepsSingleParticle>::iterator it = face_particle_impact_data_of_all_time_steps_by_threads[i].begin(); it != face_particle_impact_data_of_all_time_steps_by_threads[i].end(); ++it) {
            if (mFaceParticleImpactDataOfAllTimeSteps.find(it->first) == mFaceParticleImpactDataOfAllTimeSteps.end() ) {
                mFaceParticleImpactDataOfAllTimeSteps[it->first] = it->second;
            } 
        }
    }

    mInterParticleImpactDataOfAllParticles.push_back(time_step_database);
    mFaceParticleImpactDataOfAllParticles.push_back(face_time_step_database);

    auto d1 = GetParticleDataBase(1, inter_particle_impact_data_of_all_time_steps_by_threads[0]);
    auto dumb1 = GetParticleDataBase(1, mInterParticleImpactDataOfAllTimeSteps);
    auto dumb2 = GetParticleDataBase(2, mInterParticleImpactDataOfAllTimeSteps);    

}

void AnalyticParticleWatcher::SetNodalMaxImpactVelocities(ModelPart& analytic_model_part)
{

    #pragma omp parallel for if(analytic_model_part.NumberOfElements(0)>100)
    for (int k=0; k<(int)analytic_model_part.NumberOfElements(0); k++) {
        ElementsIteratorType i_elem = analytic_model_part.ElementsBegin() + k;
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        double& current_max_normal_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(NORMAL_IMPACT_VELOCITY);
        double& current_max_tangential_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(TANGENTIAL_IMPACT_VELOCITY);
        
        const int id = int(i_elem->Id());
        InterParticleImpactDataOfAllTimeStepsSingleParticle& particle_database = GetParticleDataBase(id, mInterParticleImpactDataOfAllTimeSteps);
        double db_normal_impact_velocity = 0.0;
        double db_tangential_impact_velocity = 0.0;
        particle_database.GetMaxCollidingSpeedFromDatabase(db_normal_impact_velocity, db_tangential_impact_velocity);
        // choose max between current and database
        current_max_normal_velocity = db_normal_impact_velocity;
        current_max_tangential_velocity = db_tangential_impact_velocity;
    }
}

void AnalyticParticleWatcher::SetNodalMaxFaceImpactVelocities(ModelPart& analytic_model_part)
{
    #pragma omp parallel for if(analytic_model_part.NumberOfElements(0)>100)
    for (int k=0; k<(int)analytic_model_part.NumberOfElements(0); k++) {
        ElementsIteratorType i_elem = analytic_model_part.ElementsBegin() + k;
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        double& current_max_normal_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(FACE_NORMAL_IMPACT_VELOCITY);
        double& current_max_tangential_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(FACE_TANGENTIAL_IMPACT_VELOCITY);

        const int n_collisions = particle.GetNumberOfCollisionsWithFaces();
        if (n_collisions){
            const int id = int(i_elem->Id());
            FaceParticleImpactDataOfAllTimeStepsSingleParticle& particle_database = GetParticleFaceDataBase(id, mFaceParticleImpactDataOfAllTimeSteps);
            double db_normal_impact_velocity = 0.0;
            double db_tangential_impact_velocity = 0.0;
            particle_database.GetMaxCollidingSpeedFromDatabase(db_normal_impact_velocity, db_tangential_impact_velocity);

            // choose max between current and database
            current_max_normal_velocity = std::max(current_max_normal_velocity, db_normal_impact_velocity);
            current_max_tangential_velocity = std::max(current_max_tangential_velocity, db_tangential_impact_velocity);
        } 
    }
}


void AnalyticParticleWatcher::SetNodalMaxLinearImpulse(ModelPart& analytic_model_part)
{
    #pragma omp parallel for if(analytic_model_part.NumberOfElements(0)>100)
    for (int k=0; k<(int)analytic_model_part.NumberOfElements(0); k++) {
        ElementsIteratorType i_elem = analytic_model_part.ElementsBegin() + k;
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));
        
        double& current_max_linear_impulse = particle.GetGeometry()[0].FastGetSolutionStepValue(LINEAR_IMPULSE);

        const int n_collisions = particle.GetNumberOfCollisions();
        if (n_collisions){
            const int id = int(i_elem->Id());
            InterParticleImpactDataOfAllTimeStepsSingleParticle& particle_database = GetParticleDataBase(id, mInterParticleImpactDataOfAllTimeSteps);
            double db_linear_impulse = 0.0;
            particle_database.GetMaxLinearImpulseFromDatabase(db_linear_impulse);            

            // choose max between current and database
            current_max_linear_impulse = std::max(current_max_linear_impulse, db_linear_impulse);
        } 
    }
}

void AnalyticParticleWatcher::ClearList(boost::python::list& my_list)
{
    while(len(my_list)){
        my_list.pop(); // only way I found to remove all entries
    }
}

void AnalyticParticleWatcher::GetParticleData(int id,
                                              boost::python::list times,
                                              boost::python::list neighbour_ids,
                                              boost::python::list normal_relative_vel,
                                              boost::python::list tangential_relative_vel)
{
    mInterParticleImpactDataOfAllTimeSteps[id].FillUpPythonLists(times, neighbour_ids, normal_relative_vel, tangential_relative_vel);
}

void AnalyticParticleWatcher::GetAllParticlesData(ModelPart& analytic_model_part,
                                                  boost::python::list times,
                                                  boost::python::list neighbour_ids,
                                                  boost::python::list normal_relative_vel,
                                                  boost::python::list tangential_relative_vel)
{
    ClearList(times);
    ClearList(neighbour_ids);
    ClearList(normal_relative_vel);
    ClearList(tangential_relative_vel);

    for (ElementsIteratorType i_elem = analytic_model_part.ElementsBegin(); i_elem != analytic_model_part.ElementsEnd(); ++i_elem){
        boost::python::list times_i;
        boost::python::list neighbour_ids_i;
        boost::python::list normal_relative_vel_i;
        boost::python::list tangential_relative_vel_i;
        const int id = int(i_elem->Id());
        GetParticleData(id, times_i, neighbour_ids_i, normal_relative_vel_i, tangential_relative_vel_i);
        times.append(times_i);
        neighbour_ids.append(neighbour_ids_i);
        normal_relative_vel.append(normal_relative_vel_i);
        tangential_relative_vel.append(tangential_relative_vel_i);
    }

}

void AnalyticParticleWatcher::GetTimeStepsData(boost::python::list ids,
                                               boost::python::list neighbour_ids,
                                               boost::python::list normal_relative_vel,
                                               boost::python::list tangential_relative_vel)
{
    ClearList(ids);
    ClearList(neighbour_ids);
    ClearList(normal_relative_vel);
    ClearList(tangential_relative_vel);
    const int n_time_steps = mInterParticleImpactDataOfAllParticles.size();

    for (int i = 0; i < n_time_steps; ++i){
        boost::python::list ids_i;
        boost::python::list neighbour_ids_i;
        boost::python::list normal_relative_vel_i;
        boost::python::list tangential_relative_vel_i;
        mInterParticleImpactDataOfAllParticles[i].FillUpPythonLists(ids_i, neighbour_ids_i, normal_relative_vel_i, tangential_relative_vel_i);
        ids.append(ids_i);
        neighbour_ids.append(neighbour_ids_i);
        normal_relative_vel.append(normal_relative_vel_i);
        tangential_relative_vel.append(tangential_relative_vel_i);
    }
}


AnalyticParticleWatcher::InterParticleImpactDataOfAllTimeStepsSingleParticle& AnalyticParticleWatcher::GetParticleDataBase(int id, std::map<int, InterParticleImpactDataOfAllTimeStepsSingleParticle>& data_base)
{
    if (mInterParticleImpactDataOfAllTimeSteps.find(id) == mInterParticleImpactDataOfAllTimeSteps.end()){
        #ifdef KRATOS_DEBUG
        KRATOS_ERROR_IF(&data_base == &mInterParticleImpactDataOfAllTimeSteps)<< "Adding an element to the database is not safe here. Should have been done before!! "<< std::endl;
        #endif 
        AnalyticParticleWatcher::InterParticleImpactDataOfAllTimeStepsSingleParticle new_particle_database(id);
        data_base[id] = new_particle_database;
        return data_base[id];
    }
    else{
        return mInterParticleImpactDataOfAllTimeSteps[id];
    }
}

AnalyticParticleWatcher::FaceParticleImpactDataOfAllTimeStepsSingleParticle& AnalyticParticleWatcher::GetParticleFaceDataBase(int id, std::map<int, FaceParticleImpactDataOfAllTimeStepsSingleParticle>& data_base)
{
    if (mFaceParticleImpactDataOfAllTimeSteps.find(id) == mFaceParticleImpactDataOfAllTimeSteps.end()){
        #ifdef KRATOS_DEBUG
        KRATOS_ERROR_IF(&data_base == &mFaceParticleImpactDataOfAllTimeSteps)<< "Adding an element to the database is not safe here. Should have been done before!! "<< std::endl;
        #endif 
        AnalyticParticleWatcher::FaceParticleImpactDataOfAllTimeStepsSingleParticle new_particle_database(id);
        data_base[id] = new_particle_database;
        return data_base[id];
    }
    else{
        return mFaceParticleImpactDataOfAllTimeSteps[id];
    }
}

/// Turn back information as a string.
std::string AnalyticParticleWatcher::Info() const {
        return "AnalyticParticleWatcher";
}

/// Print information about this object.
void AnalyticParticleWatcher::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AnalyticParticleWatcher::PrintData(std::ostream& rOStream) const {}


} // namespace Kratos
