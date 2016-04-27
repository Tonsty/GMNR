#ifndef GMNR_TIME_H
#define GMNR_TIME_H

#ifdef WIN32

# include <limits.h>
# include <windows.h>

namespace gmnr {

  struct timestamp { LARGE_INTEGER t; };

  static inline double LI2d(const LARGE_INTEGER &li)
  {
	// Work around random compiler bugs...
	double d = *(unsigned *)(&(li.HighPart));
	d *= 65536.0 * 65536.0;
	d += *(unsigned *)(&(li.LowPart));
	return d;
  }

  static inline float operator - (const timestamp &t1, const timestamp &t2)
  {
	static LARGE_INTEGER PerformanceFrequency;
	static int status = QueryPerformanceFrequency(&PerformanceFrequency);
	if (status == 0) return 1.0f;

	return (LI2d(t1.t) - LI2d(t2.t)) / LI2d(PerformanceFrequency);
  }

  static inline timestamp now()
  {
	timestamp t;
	QueryPerformanceCounter(&t.t);
	return t;
  }

}

#else

# include <sys/time.h>
# include <unistd.h>

namespace gmnr {

  typedef struct timeval timestamp;

  static inline float operator - (const timestamp &t1, const timestamp &t2)
  {
	return (float)(t1.tv_sec  - t2.tv_sec) +
	       1.0e-6f*(t1.tv_usec - t2.tv_usec);
  }

  static inline timestamp now()
  {
	timestamp t;
	gettimeofday(&t, 0);
	return t;
  }

}

#endif

#endif
