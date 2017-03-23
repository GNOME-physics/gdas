Extracting GPS time range
=========================

We use the `LIGOTimeGPS <http://software.ligo.org/docs/lalsuite/lal/group___l_a_l_datatypes.html#ss_LIGOTimeGPS>`_ structure from the =glue.lal= package to /store the starting and ending time in the dataset to nanosecond precision and synchronized to the Global Positioning System time reference/. Once both times are defined, the range of value is stored in a semi-open interval using the `segment <http://software.ligo.org/docs/glue/glue.__segments.segment-class.html>`_ module from the =glue.segments= package. ::

  # Starting epoch relative to GPS starting epoch
  start_time = LIGOTimeGPS(args.analysis_start_time or args.gps_start_time)
  # Ending epoch relative to GPS ending epoch
  end_time = LIGOTimeGPS(args.analysis_end_time or args.gps_end_time)
  # Represent the range of values in the semi-open interval 
  inseg = segment(start_time,end_time)

Prepare output file for given time range
========================================

::

  xmldoc = ligolw.Document()
  xmldoc.appendChild(ligolw.LIGO_LW())
  
  ifo = args.channel_name.split(":")[0]
  proc_row = register_to_xmldoc(xmldoc, __program__, args.__dict__, ifos=[ifo],version=glue.git_version.id, cvs_repository=glue.git_version.branch, cvs_entry_time=glue.git_version.date)
  
  # Figure out the data we actually analyzed
  outseg = determine_output_segment(inseg, args.psd_segment_length, args.sample_rate, window_fraction)
  
  ss = append_search_summary(xmldoc, proc_row, ifos=(station,), inseg=inseg, outseg=outseg)
  
  for sb in event_list:
      sb.process_id = proc_row.process_id
      sb.search = proc_row.program
      #sb.ifo, sb.channel = args.channel_name.split(":")
      sb.ifo, sb.channel = station, setname
  
  xmldoc.childNodes[0].appendChild(event_list)
  fname = make_filename(station, inseg)
  
  utils.write_filename(xmldoc, fname, gz=fname.endswith("gz"), verbose=True)

Plot trigger results
====================

::

  events = SnglBurstTable.read(fname+'.gz')
  #del events[10000:]
  plot = events.plot('time', 'central_freq', "duration", "bandwidth", color="snr")
  #plot = events.plot('time', 'central_freq', color='snr')
  #plot.set_yscale("log")
  plot.set_ylim(1e-0, 250)
  t0 = 1153742417
  plot.set_xlim(t0 + 0*60, t0 + 1*60)
  #plot.set_xlim(t0 + 28, t0 + 32)
  pyplot.axvline(t0 + 30, color='r')
  cb = plot.add_colorbar(cmap='viridis')
  plot.savefig("triggers.png")
