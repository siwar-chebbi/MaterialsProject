from pymatgen import MPRester

mpr = MPRester("78OAi0lR9kdkyiAi")


mpids = ['mp-1021516', 'mp-9580'] # replace with your list of missing material ids

results = mpr.query({"task_ids": {"$in": mpids}}, ["task_id"])

assert len(results) == len(mpids)