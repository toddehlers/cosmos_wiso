#if defined use_comm_SIPC || defined use_comm_GMEM
#include <signal.h>
extern void getfpe_();
void fpecatch_()
{   signal(SIGFPE,getfpe_); }
#endif
