#ifdef use_comm_SIPC
#include <signal.h>
extern void ferror_();
void cldcatch_()
{   signal(SIGCLD,ferror_); }
#endif
