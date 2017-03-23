.. _tileduration:

Explore multiple tile durations
===============================

Now that we create a tile with a specific bandwidth, we can start exploring different durations for the tile. We will start checking if the user manually defined a value for the longest duration tile to compute, which can be done using the ``--max-duration`` argument. If not, the value will be set to 32. ::

          if args.max_duration is not None:
              max_dof = 2 * args.max_duration * (band * (nc_sum+1))
          else:
              max_dof = 32
          assert max_dof >= 2

Since we produce (initially) tiles with 1 degree of freedom, the duration goes as one over twice the bandwidth. ::

          print "\t\t...getting longer durations..."
          #for j in [2**l for l in xrange(1, int(math.log(max_dof, 2))+1)]:
          for j in [2**l for l in xrange(0, int(math.log(max_dof, 2)))]:
              sys.stderr.write("\t\tSumming DOF = %d ..." % (2*j))
              #tlen = tiles.shape[1] - j + 1
              tlen = tiles.shape[1] - 2*j + 1 + 1
              if tlen <= 0:
                  print >>sys.stderr, " ...not enough samples."
                  continue
              dof_tiles = numpy.zeros((tiles.shape[0], tlen))
              #:sum_filter = numpy.ones(j)
              # FIXME: This is the correct filter for 50% overlap
              sum_filter = numpy.array([1,0] * (j-1) + [1])
              #sum_filter = numpy.array([1,0] * int(math.log(j, 2)-1) + [1])
              for f in range(tiles.shape[0]):
                  # Sum and drop correlate tiles
                  # FIXME: don't drop correlated tiles
                  #output = numpy.convolve(tiles[f,:], sum_filter, 'valid')
                  dof_tiles[f] = fftconvolve(tiles[f], sum_filter, 'valid')
              print >>sys.stderr, " done"
              print "Summed tile energy mean: %f, var %f" % (numpy.mean(dof_tiles), numpy.var(dof_tiles))
              level_tdiff = time.time() - tdiff
              print >>sys.stderr, "Done with this resolution, total %f" % level_tdiff

Finally, the bandwidth and duration of the tile can be defined as followed: ::

              # Current bandwidth of the time-frequency map tiles
              current_band = band * (nc_sum + 1)
              # How much each "step" is in the frequency domain -- almost
              # assuredly the fundamental bandwidth
              df = current_band
              # How much each "step" is in the time domain -- under sampling rate
              # FIXME: THis won't work if the sample rate isn't a power of 2
              dt = 1.0 / 2 / (2 * current_band) * 2
              full_band = 250
              dt = current_band / full_band * ts_data.sample_rate
              dt = 1.0/dt
              # Duration is fixed by the NDOF and bandwidth
              duration = j / 2.0 / current_band    
