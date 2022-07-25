////////////////////////////////////////////////////////////////////////
// Class:       AbsRunningSumTPFinderPass1
// File:        AbsRunningSumTPFinderPass1_service.cc
////////////////////////////////////////////////////////////////////////

#include "duneana/DAQSimAna/AbsRunningSumHitFinder/AbsRunningSumTPFinderPass1.h"
#include "duneana/DAQSimAna/AbsRunningSumHitFinder/AbsRSAlgParts.h"
#include "duneana/DAQSimAna/AlgParts.h"

#include <algorithm> // for std::transform
#include <numeric> // for std::accumulate


AbsRunningSumTPFinderPass1::AbsRunningSumTPFinderPass1(fhicl::ParameterSet const & p)
  : m_threshold          (p.get<unsigned int>      ("Threshold"            )),
    m_useSignalKill      (p.get<bool>              ("UseSignalKill"        ,                        true)),
    m_signalKillLookahead(p.get<short>             ("SignalKillLookahead"  ,                           5)),
    m_signalKillThreshold(p.get<short>             ("SignalKillThreshold"  ,                          15)),
    m_signalKillNContig  (p.get<short>             ("SignalKillNContig"    ,                           1)),
    m_frugalNContig      (p.get<short>             ("FrugalPedestalNContig",                          10)),
    m_doFiltering        (p.get<bool>              ("DoFiltering"          ,                        true)),
    m_downsampleFactor   (p.get<unsigned int>      ("DownsampleFactor"     ,                           1)),
    m_filterTaps         (p.get<std::vector<short>>("FilterCoeffs"         , {2,  9, 23, 31, 23,  9,  2})),
    m_multiplier         (std::accumulate(m_filterTaps.begin(), m_filterTaps.end(), 0)),
    m_R                  (p.get<float>             ("R"      ,                     0.7))
   
{

}

std::vector<short> AbsRunningSumTPFinderPass1::downSample(const std::vector<short>& orig) {

  //---------------------------------------------
  // Do the downsampling
  //---------------------------------------------
  if (m_downsampleFactor==1) {
    return orig;
  }
  else {
    std::vector<short> waveform;
    for(size_t i=0; i<orig.size(); i+=m_downsampleFactor) {
      waveform.push_back(orig[i]);
    }
    return waveform;
  }
}

std::vector<short> AbsRunningSumTPFinderPass1::findPedestal(const std::vector<short>& waveform)
{
  //---------------------------------------------
  // Pedestal subtraction
  //---------------------------------------------
  const std::vector<short>& pedestal=m_useSignalKill ?
    frugal_pedestal_sigkill(waveform,
			    m_signalKillLookahead,
			    m_signalKillThreshold,
			    m_signalKillNContig) :
    frugal_pedestal(waveform, m_frugalNContig);
  return pedestal;
}

std::vector<short> AbsRunningSumTPFinderPass1::filter(const std::vector<short>& pedsub) {

  //---------------------------------------------
  // Filtering
  //---------------------------------------------
  const size_t ntaps = m_filterTaps.size();
  const short*  taps = m_filterTaps.data();

  std::vector<short> filtered(m_doFiltering ? 
			      apply_fir_filter(pedsub, ntaps, taps) :
			      pedsub);
  if (!m_doFiltering) {
    std::transform(filtered.begin(), filtered.end(),
		   filtered.begin(), 
		   [=](short a) { return a*m_multiplier; });
  }
  return filtered;
}

std::vector<short> AbsRunningSumTPFinderPass1::AbsRunningSum(const std::vector<short>& filtered, float R) {

  //---------------------------------------------
  // Absolute Running Sum Algorithm
  //---------------------------------------------

  //initialise 
  short s = 2;
  std::vector<short> absRS(filtered.size(), 0); absRS[0] = filtered[0]/s;
  for (size_t i=0; i<filtered.size(); ++i) { 
    absRS[i] = R*absRS[i-1] + std::abs(filtered[i]/s);
  }
  return absRS;
}


std::vector<short> AbsRunningSumTPFinderPass1::findPedestal_absRS(const std::vector<short>& waveform) {

  //---------------------------------------------
  // absRS pedestal subtraction
  //---------------------------------------------

  const std::vector<short>& pedestal = m_useSignalKill ? 
    frugal_pedestal_sigkill(waveform, m_signalKillLookahead, m_signalKillThreshold, m_signalKillNContig) : 
    frugal_pedestal_absRS(waveform, m_frugalNContig);
  return pedestal;
}

void
AbsRunningSumTPFinderPass1::hitFinding(const std::vector<short>& waveform,
				    std::vector<AbsRunningSumTPFinderTool::Hit>& hits,
				    int channel) {

  //---------------------------------------------
  // Hit finding
  //---------------------------------------------
  bool is_hit  = false;
  bool was_hit = false;
  AbsRunningSumTPFinderTool::Hit hit(channel, 0, 0, 0);
  for(size_t isample=0; isample<waveform.size()-1; ++isample){
    int   sample_time = isample * m_downsampleFactor;
    short adc         = waveform[isample];
    is_hit = adc >  (short)m_threshold;
    if(is_hit && !was_hit) {
      hit.startTime         = sample_time;
      hit.charge            = adc;
      hit.timeOverThreshold = m_downsampleFactor;
    }
    if(is_hit && was_hit) {
      hit.charge            += adc*m_downsampleFactor;
      hit.timeOverThreshold += m_downsampleFactor;
    }
    if(!is_hit && was_hit) {
      hit.charge /= m_multiplier;
      hits.push_back(hit);
    }
    was_hit = is_hit;
  }
}

std::vector<AbsRunningSumTPFinderTool::Hit>
AbsRunningSumTPFinderPass1::findHits(const std::vector<unsigned int>& channel_numbers, 
				  const std::vector<std::vector<short>>& adc_samples) {

  auto hits = std::vector<AbsRunningSumTPFinderTool::Hit>();
  std::cout << "findHits called with "      << adc_samples.size()
	    << " channels. First chan has " << adc_samples[0].size() << " samples" << std::endl;

  for(size_t ich=0; ich<adc_samples.size(); ++ich){
    const std::vector<short>& waveformOrig = adc_samples[ich];

    std::vector<short> waveform  = downSample  (waveformOrig);
    std::vector<short> pedestal  = findPedestal(waveform    );
    std::vector<short> pedsub(waveform.size(), 0);
    for(size_t i=0; i<pedsub.size(); ++i)
      pedsub[i]=waveform[i]-pedestal[i];
    std::vector<short> filtered  = filter(pedsub);
    std::vector<short> AbsRS = AbsRunningSum(filtered, m_R);
    std::vector<short> AbsRS_ped = findPedestal_absRS(AbsRS);
    std::vector<short> AbsRS_pedsub(AbsRS.size(), 0); 
    for (size_t j = 0; j < AbsRS_pedsub.size(); ++j){
      AbsRS_pedsub[j] = AbsRS[j] - AbsRS_ped[j];
    }
  
    hitFinding(AbsRS_pedsub, hits, channel_numbers[ich]);
  }
  std::cout << "Returning " << hits.size() << " hits for a threshold of " << m_threshold << std::endl;
  return hits;
}

DEFINE_ART_CLASS_TOOL(AbsRunningSumTPFinderPass1)
