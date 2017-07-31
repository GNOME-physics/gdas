def coincidence():
    # 1. Choose a chunk of time that will complete in a reasonable amount of time.
    gps_start_time=XXXX
    gps_end_time=YYYY
    ## Loop on instruments:
    #for inst in "B1", "K1", "F1" ...; do
    #    # 2. Run pyburst_excesspower_gnome on that chunk of time on that instrument
    #     pyburst_excesspower_gnome ...
    #    # 3. Run lalapps_bucluster on the output from excesspower.
    #    lalapps_bucluster -v -p pyburst_excesspower_gnome -c excesspower (input files)
    #done
    ## 4. Add together all the results from 2 & 3
    #ligolw_add -v (clustered input files) -o ouput_name.xml.gz
    ## 5. Run lalapps_burca on the output from bucluster.
    #lalapps_burca -v output_name.xml.gz
