from __future__ import division
from libtbx import easy_run
import sys, os

num_procs = 1
def do_work(path):
  root        = os.path.dirname(path)
  basename    = os.path.splitext(os.path.basename(path))[0]
  datablock   = os.path.join(root, basename + "_datablock.json")
  spots       = os.path.join(root, basename + "_strong.pickle")
  reflections = os.path.join(root, basename + "_reflections.pickle")
  experiments = os.path.join(root, basename + "_experiments.json")

  print "Preparing to index", basename

  easy_run.call("dials.import %s --output %s"%(path,datablock))
  easy_run.call("dials.find_spots %s threshold.xds.sigma_strong=15 min_spot_size=2 -o %s"%(datablock, spots))
  easy_run.call("dials.index %s %s method=fft1d beam.fix=all detector.fix=all known_symmetry.unit_cell=93,93,130,90,90,120 known_symmetry.space_group=P6122 n_macro_cycles=5 d_min_final=0.5 experiments=%s reflections=%s"%(spots, datablock, experiments, reflections))

  if os.path.exists(reflections):
    assert os.path.exists(experiments)
    print basename, "indexed succesfully"

    results.append((reflections,experiments))
  else:
    print basename, "failed to index"

if num_procs > 1:
  from multiprocessing import Manager, Process

  mgr = Manager()
  results = mgr.list()
  def worker():
    for item in iter( q.get, None ):
      do_work(item)
      q.task_done()
    q.task_done()

  q = mgr.JoinableQueue()
  procs = []
  for i in range(num_procs):
    procs.append( Process(target=worker) )
    procs[-1].daemon = True
    procs[-1].start()

  for path in sys.argv[1:]:
    q.put(path)

  q.join()

  for p in procs:
    q.put( None )

  q.join()

  for p in procs:
    p.join()
else:
  results = []
  for path in sys.argv[1:]:
    do_work(path)


print "Writing phil input to metrology refiner..."

f = open("all.phil", 'w')
for (indexed, experiments) in results:
  f.write("input {\n")
  f.write("  experiments = %s\n"%experiments)
  f.write("  reflections = %s\n"%indexed)
  f.write("}\n")
f.close()

print "Done"

