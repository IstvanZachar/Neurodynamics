#include <stdio.h>
#include <windows.h>

// Some ANSI escape sequences.
// These only work if the terminal the code is run on supports ANSI escapes (mostly Unix/Linux systems, but not Windows
// for more escape sequences, refer to: http://ascii-table.com/ansi-escape-sequences.php


#define STYLE_BLACK    "\x1B[30m"
#define STYLE_RED      "\x1B[31m"
#define STYLE_GREEN    "\x1B[32m"
#define STYLE_YELLOW   "\x1B[33m"
#define STYLE_BLUE     "\x1B[34m"
#define STYLE_MAGENTA  "\x1B[35m"
#define STYLE_CYAN     "\x1B[36m"
#define STYLE_WHITE    "\x1B[37m"
#define STYLE_FG_RESET "\x1B[0m"

#define STYLE_BG_RED   "\033[41m"
#define STYLE_BG_RESET "\033[49m"

#define STYLE_RESET    "\033[39;49m" // reset FG and BG
#define STYLE_BOLD     "\x1B[1m"

// USAGE
// printf(STYLE_RED "This text is red." STYLE_RESET "\n");



// If terminal does not support ANSI escape codes, access the OS's console directly, via "windows.h" functions
// I've only implemented red and green printing functions.
// These should be used the same way as printf, with any number and types of arguments
// style integers:
//   -1: reset
//   (1, 15) black BG, colored FG
//   (16, 255) same FG colors as above, with different BG colors

void printRed(char* fmt, ...) {
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12);
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // reset
}

void printGreen(char* fmt, ...) {
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 2);
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // reset
}

void printCyan(char* fmt, ...) {
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 11);
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // reset
}

void printStyle(int style, char* fmt, ...) {
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), style);
	va_list args;
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7);
}

