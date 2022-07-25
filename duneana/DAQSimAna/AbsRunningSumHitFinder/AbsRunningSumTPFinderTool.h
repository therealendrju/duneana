#ifndef AbsRunningSumTPFinderTool_h
#define AbsRunningSumTPFinderTool_h

#include <vector>
#include <iostream>

class AbsRunningSumTPFinderTool {
 
 public:
  struct Hit
  {
  Hit(int _channel, int _startTime, int _charge, int _timeOverThreshold)
  : channel(_channel),
      startTime(_startTime),
      charge(_charge),
      timeOverThreshold(_timeOverThreshold)
    {}
    int channel;
    int startTime;
    int charge;
    int timeOverThreshold;
  };

  virtual ~AbsRunningSumTPFinderTool() =default;

  virtual std::vector<AbsRunningSumTPFinderTool::Hit>
    findHits(const std::vector<unsigned int>& channel_numbers, 
             const std::vector<std::vector<short>>& adc_samples) = 0;
 
};

#endif // include guard
