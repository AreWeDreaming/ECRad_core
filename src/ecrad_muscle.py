from libmuscle import Instance, Message
from ymmsl import Operator
import logging
import imas
import os
import yaml
import time
import sys

def ECRad_MUSCLE3_test():
    """
    Regression test for MUSLCE3 ECRad
    """
    logging.info("ECRad MUSCLE3 test start!")
    instance = Instance({
            Operator.O_I:    ['ECRad_report'],
            Operator.S:      ['ECRad_task']})
    while instance.reuse_instance():
        scenario_path = instance.get_setting('scenario_info', 'str')
        scenario = os.path.join(os.path.dirname(__file__), scenario_path)
        with open(scenario, 'r') as scenario_file:
            config = yaml.load(scenario_file,Loader=yaml.CLoader)
        ids = {}
        for ids_id in ["equilibrium", "core_profiles", "ece", "wall"]:
            logging.info(f"Loading: {ids_id} ids")
            input = imas.DBEntry(imas.imasdef.MDSPLUS_BACKEND, config['db_' + ids_id],
                              config['shot_' + ids_id], config['run_' + ids_id], 
                              config['user_' + ids_id])
            input.open()
            if ids_id in ["equilibrium", "core_profiles"]:
                ids[ids_id] = input.get_slice(ids_id,config['time_slice_'+ ids_id],1)
            else:
                ids[ids_id] = input.get(ids_id)
            input.close()
        logging.info("Test sending INIT task")

        msg = Message(time.time(), data=["INIT", ids["wall"].serialize(), ids["equilibrium"].serialize(), ids["ece"].serialize()])
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        logging.info("Test received report from ECRad")
        if msg.data[0] != "Init success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        logging.info("Test sending timepoint task")
        msg = Message(time.time(), data=["Timepoint", 1])
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        logging.info("Test received report from ECRad")
        if msg.data[0] != "Timepoint success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        msg = Message(time.time(), data=["Run", ids["core_profiles"]])
        logging.info("Test sending run task")
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        logging.info("Sending run task")
        if msg.data[0] != "Run success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        break
    

if __name__ == '__main__':
    logging.basicConfig(filename='ecrad_tester.log', level=logging.INFO)
    logging.getLogger().setLevel(logging.INFO)
    try:
        ECRad_MUSCLE3_test()
    except Exception as e:
        logging.error(str(e))
    

    

