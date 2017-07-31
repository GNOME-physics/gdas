import math
from glue.segments import segment
from glue.lal import LIGOTimeGPS
from glue.ligolw import lsctables, utils, ligolw
lsctables.use_in(ligolw.LIGOLWContentHandler)
from glue.ligolw.utils import process, search_summary
from scipy.stats import poisson, uniform, expon, chi2

def fake_trigger_generator(instrument='H1'):

    xmldoc = ligolw.Document()
    xmldoc.appendChild(ligolw.LIGO_LW())
    
    # Process information
    proc = process.append_process(xmldoc, "fake_search")
    process.append_process_params(xmldoc, proc, {})
    
    t0 = 1e9
    ntrig = 1000
    ifo = instrument
    inseg = segment(LIGOTimeGPS(t0), LIGOTimeGPS(t0 + ntrig / 10))
    outseg = segment(LIGOTimeGPS(t0), LIGOTimeGPS(t0 + ntrig / 10))
    
    # Search summary
    search_summary.append_search_summary(xmldoc, proc, comment="Fake triggers", ifos=(ifo,), inseg=inseg, outseg=outseg)
    
    columns = ['chisq_dof', 'bandwidth', 'central_freq', 'confidence', 'peak_time_ns', 'start_time',
               'process_id', 'fhigh', 'stop_time_ns', 'channel', 'ifo', 'duration', 'event_id', 'hrss',
               'stop_time', 'peak_time', 'snr', 'search', 'start_time_ns', 'flow', 'amplitude']
    table = lsctables.New(lsctables.SnglBurstTable, columns)
    
    # Generate uniformly distributed trigger times with approximate rate of 10 s
    times = t0 + uniform.rvs(0, ntrig / 10., ntrig)
    
    for t in times:
        row = table.RowType()
    
        # time frequency position and extent
        row.chisq_dof = int(2 + expon.rvs(2))
        row.duration = 1. / 2**int(uniform.rvs(0, 7))
        row.bandwidth = row.chisq_dof / row.duration / 2
    
        row.central_freq = uniform.rvs(16, 2048)
        row.flow = max(row.central_freq - row.bandwidth, 0)
        row.fhigh = min(row.central_freq + row.bandwidth, 2048)
    
        ns, sec = math.modf(t)
        ns = int("%09d" % (ns * 1e9))
        row.peak_time, row.peak_time_ns = int(sec), ns
    
        ns, sec = math.modf(t - row.duration / 2)
        ns = int("%09d" % (ns * 1e9))
        row.start_time, row.start_time_ns = int(sec), ns
    
        ns, sec = math.modf(t + row.duration / 2)
        ns = int("%09d" % (ns * 1e9))
        row.stop_time, row.stop_time_ns = int(sec), ns
    
        # TODO: Correlate some triggers, an upward fluctuation often triggers a few
        # tiles ontop of each other
    
        # SNR and confidence
        row.snr = 5.
        while row.snr < 2 * row.chisq_dof:
            row.snr = chi2.rvs(row.chisq_dof)
        row.confidence = chi2.sf(row.snr, row.chisq_dof)
        row.snr = math.sqrt(row.snr / row.chisq_dof - 1)
        row.hrss = row.amplitude = 1e-21
    
        # metadata
        row.search = "fake_search"
        row.channel = "FAKE"
        row.ifo = ifo
    
        row.event_id = table.get_next_id()
        row.process_id = proc.process_id
    
        table.append(row)
    
    xmldoc.childNodes[0].appendChild(table)
    
    utils.write_filename(xmldoc, "%s-FAKE_SEARCH-%d-%d.xml.gz" % (ifo, int(t0), 10000), gz=True)
