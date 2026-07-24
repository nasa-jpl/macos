#include <signal.h>
#include <unistd.h>

static void sigint_handler(int sig)
{
    const char msg[] = "\nMACOS: interrupted\n";
    write(STDOUT_FILENO, msg, sizeof(msg) - 1);
    _exit(0);
}

void install_sigint_handler(void)
{
    struct sigaction sa;
    sa.sa_handler = sigint_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGINT, &sa, NULL);
}
