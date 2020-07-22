from pymatgen import MPRester

mpr = MPRester("eDCEK5m9WVjmajp7e8af")

mpids = ['mp-1021516', 'mp-9580'] # replace with your list of missing material ids

results = mpr.query({"task_ids": {"$in": mpids}}, ["task_id"])
print (results)
assert len(results) == len(mpids)