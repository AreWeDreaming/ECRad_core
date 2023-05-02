from libmuscle import Instance, Message
from ymmsl import Operator
import logging
import imas
import os
import yaml
import time

def ECRad_MUSCLE3_test():
    """
    Regression test for MUSLCE3 ECRad
    """

    instance = Instance({
            Operator.F_INIT: ['ECRad_init'],
            Operator.O_I:    ['ECRad_task'],
            Operator.S:      ['ECRad_report']})
    
    while instance.reuse_instance():
        scenario_path = instance.get_setting('scenario_info', 'str')
        scenario = os.path.join(os.path.dirname(__file__), scenario_path)
        with open(scenario, 'r') as scenario_file:
            config = yaml.load(scenario_file,Loader=yaml.CLoader)
        ids = {}
        for ids_id in ["equilibrium", "core_profiles", "ece", "wall"]:
            print(f"Loading: {ids_id} ids")
            input = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, config['db_' + ids_id],
                              config['shot_' + ids_id], config['run_' + ids_id], 
                              config['user_' + ids_id])
            input.open()
            if ids_id in ["equilibrium", "core_profiles"]:
                ids[ids_id] = input.get_slice(ids_id,config['time_slice_'+ ids_id],1)
            else:
                ids[ids_id] = input.get(ids_id)
            input.close()
        print("Test sending INIT task")
        msg = Message(time.time(), data=["INIT", ids["wall"], ids["equilibrium"], ids["ece"]])
        instance.send("ECRad_init", msg)
        msg = instance.receive("ECRad_report")
        print("Test received report from ECRad")
        if msg.data[0] != "Init success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        print("Test sending timepoint task")
        msg = Message(time.time(), data=["Timepoint", 1])
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        print("Test received report from ECRad")
        if msg.data[0] != "Timepoint success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        msg = Message(time.time(), data=["Run", ids["core_profiles"]])
        print("Test sending run task")
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        print("Sending run task")
        if msg.data[0] != "Run success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        break
    

if __name__ == '__main__':
    logging.basicConfig()
    logging.getLogger().setLevel(logging.INFO)
    ECRad_MUSCLE3_test()
    
    

    

