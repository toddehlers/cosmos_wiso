#ifdef use_comm_SIPC
#include <signal.h>
extern void fignore_();
void cldignore_()
{   signal(SIGCLD,fignore_); }
#endif
