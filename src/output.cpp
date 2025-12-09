#include "output.h"

void Output::ShowWarning(const char* message)
{
    printf("\n################  WARNING ################\n");
    printf("%s\n", message);
    printf("########################################\n");
}

void Output::ShowError(const char* message)
{
    printf("\n################    ERROR   ############\n");
    printf("%s\n", message);
    printf("########################################\n");
}

void Output::ShowMessage(const char* message)
{
    printf("\n######################################\n");
    printf("%s\n", message);
    printf("########################################\n");
}

