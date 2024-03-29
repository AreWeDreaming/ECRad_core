
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
            Operator.S:    ['ECRad_report'],
            Operator.O_I:  ['ECRad_task']})
    while instance.reuse_instance():
        scenario_path = instance.get_setting('scenario_info', 'str')
        scenario = os.path.join(os.path.dirname(__file__), scenario_path)
        with open(scenario, 'r') as scenario_file:
            config = yaml.load(scenario_file,Loader=yaml.CLoader)
        ids = {}
        for ids_id in ["equilibrium", "core_profiles", "ece", "wall"]:
            logging.info(f"Loading: {ids_id} ids")
            input = imas.DBEntry(imas.imasdef.HDF5_BACKEND, config['db_' + ids_id],
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
        msg = Message(time.time(), data=["Timepoint", 1, ids["core_profiles"].serialize(), 1])
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        logging.info("Test received report from ECRad")
        if msg.data[0] != "Timepoint success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        msg = Message(time.time(), data=["Run", ids["core_profiles"].serialize()])
        logging.info("Test sending run task")
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        logging.info("Sending run task")
        if msg.data[0] != "Run success":
            raise ValueError(f"ECRad reports: {msg.data[0]}")
        logging.info("Run finished")

        logging.info("Saving results")
        output = imas.DBEntry(imas.imasdef.HDF5_BACKEND, 'ITER', config['shot_equilibrium'], config['run_out'], config['output_user'])
        output.create()
        ece_out = imas.ece()
        logging.info("Deserializing ece")
        ece_out.deserialize(msg.data[1])
        output.put(ece_out)
        logging.info("Finished saving results")
        msg = Message(time.time(), data=["Exit"])
        logging.info("Test sending run task")
        instance.send("ECRad_task", msg)
        msg = instance.receive("ECRad_report")
        logging.info("Sending run task")
        break
    

if __name__ == '__main__':
    logging.basicConfig(filename='ecrad_tester.log', level=logging.INFO)
    logging.getLogger().setLevel(logging.INFO)
    try:
        if "debug=True" in sys.argv:
            import debugpy
            debugpy.listen(5678)
            print("Waiting for debugger attach")
            debugpy.wait_for_client()
        ECRad_MUSCLE3_test()
    except Exception as e:
        logging.error(str(e))
    

    

