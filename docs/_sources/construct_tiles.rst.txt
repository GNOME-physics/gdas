.. _tilebandwidth:

Constructing tiles of different bandwidth
=========================================

First and foremost, we define a clipping region in the data to be used to remove window corruption, this is non-zero if the ``window_fraction`` variable is set to a non-zero value. ::

      print tprint(t0,t1),"Beginning tile construction..."
      # Clip the boundaries to remove window corruption
      clip_samples = int(args.psd_segment_length * window_fraction * args.sample_rate / 2)

In order to perform a multi-resolution search, tiles of many different bandwidths and durations will be scanned. We first need to setup a loop such that the maximum number of additional channel is equal to the base 2 logarithm of the total number of channels. The number of narrow band channels to be summed (``nc_sum``) would therefore be equal to 2 to the power of the current quantity of additional channels. ::

      for nc_sum in range(0, int(math.log(nchans, 2)))[::-1]: # nc_sum additional channel adds
          nc_sum = 2**nc_sum - 1
          print tprint(t0,t1,t2),"Summing %d narrow band channels..." % (nc_sum+1)

The undersampling rate for this tile can be calculated using the channel frequency band and the number of narrow band channels to be summed such that the bandwidth of the tile is equal to ``band * (nc_sum + 1)``. ::

          us_rate = int(round(1.0 / (2 * band*(nc_sum+1) * ts_data.delta_t)))
          print >>sys.stderr, "Undersampling rate for this level: %f" % (args.sample_rate/us_rate)

"Virtual" wide bandwidth channels are constructed by summing the samples from multiple channels, and correcting for the overlap between adjacent channel filters. We then define the normalised channel at the current level and create a time frequency map for this tile using the :ref:`make_indp_tiles <make_indp_tiles>` internal function. In other word, we are constructing multiple sub-tiles for which we can determined the respective energy in the given frequency band. ::

          mu_sq = mu_sq_dict[nc_sum]
          sys.stderr.write("\t...calculating tiles...")
          if clip_samples > 0:
              tiles = make_indp_tiles(tf_map[:,clip_samples:-clip_samples:us_rate], nc_sum, mu_sq)
          else:
              tiles = make_indp_tiles(tf_map[:,::us_rate], nc_sum, mu_sq)
          sys.stderr.write(" TF-plane is %dx%s samples... " % tiles.shape)
          print >>sys.stderr, " done"
          print "Tile energy mean: %f, var %f" % (numpy.mean(tiles), numpy.var(tiles))
