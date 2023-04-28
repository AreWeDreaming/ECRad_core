from libmuscle import Instance, Grid, Message
from ymmsl import Operator
import random,copy
import logging
import imas

def fake_integrator() -> None:
    """Muscled fake integrator
    """

    instance = Instance({
            Operator.F_INIT: ['core_profiles_in'],
            Operator.O_I:    ['core_profiles_out'],
            Operator.S:      ['dynamic_interferometer_in'],
            Operator.O_F:    ['final_core_profiles_out','final_interferometer_out']})

    FirstRun = True
    
    while instance.reuse_instance():
        
        # begin F_INIT
        msg_core_profiles = instance.receive("core_profiles_in")
        core_profiles = imas.core_profiles()
        core_profiles.deserialize(msg_core_profiles.data)
        t_core_profiles = msg_core_profiles.timestamp
        # end F_INIT
        
        n_iterations = 10
        for i_iteration in range(n_iterations):
    
            # begin O_I
            core_profiles_out = copy.deepcopy(core_profiles)
            core_profiles_out.profiles_1d[0].electrons.density = \
                copy.deepcopy(core_profiles.profiles_1d[0].electrons.density*random.random())
            t_cur = core_profiles_out.time[-1]
            msg = Message(t_cur, data=core_profiles_out.serialize())
            instance.send("core_profiles_out", msg)
            # end O_I

            # begin S
            msg_dynamic_interferometer_in = instance.receive("dynamic_interferometer_in")
            dynamic_interferometer_in = imas.interferometer()
            dynamic_interferometer_in.deserialize(msg_dynamic_interferometer_in.data)
            t_dynamic_interferometer_in = msg_dynamic_interferometer_in.timestamp
            # end S

        t_cur = core_profiles_out.time[-1]
        msg = Message(t_cur)
        instance.send("core_profiles_out", msg)
                    
        # begin O_F
        t_cur = core_profiles_out.time[-1]
        msg = Message(t_cur, data=core_profiles_out.serialize())
        instance.send("final_core_profiles_out", msg)
        t_cur = dynamic_interferometer_in.time[-1]
        msg = Message(t_cur, data=dynamic_interferometer_in.serialize())
        instance.send("final_interferometer_out", msg)
        # end O_F
    
if __name__ == '__main__':
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    fake_integrator()
    

    

