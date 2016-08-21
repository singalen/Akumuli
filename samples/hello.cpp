#include <cstdlib>
#include <iostream>
#include "akumuli.h"

/** Logger f-n that shuld be used in libakumuli */
static void static_logger(aku_LogLevel tag, const char * msg) {
//	static Logger logger = Logger("Main", 32);
	switch(tag) {
		case AKU_LOG_ERROR:
			std::cerr << msg << "\n";
			break;
		case AKU_LOG_INFO:
			std::cout << msg << "\n";
			break;
		case AKU_LOG_TRACE:
			std::cout << msg << "\n";
			break;
	}
}

/** Panic handler for libakumuli.
  * Shouldn't be called directly, writes error message and
  * writes coredump (this depends on system configuration)
  */
void panic_handler(const char * msg) {
	// write error message
	static_logger(AKU_LOG_ERROR, msg);
	static_logger(AKU_LOG_ERROR, "Terminating (core dumped)");
	// this should generate SIGABORT and triger coredump
	abort();
}


void writeDummy(aku_Database* db) {
	const int count = 10000;
	std::cout << "Writing dummy data\n";
	for(uint64_t i = 0; i < count; i++) {
//		int someBuffer = 1;
//		aku_MemRange memr;
//		memr.address = (void*) &someBuffer;
//		memr.length = sizeof(someBuffer);
		
		aku_Sample sample;
		sample.timestamp = i;
		sample.paramid = 42;
//		sample.payload = memr;
		
		if (i % (count / 50) == 0) {
			std::cout << '.';
		}
		
		sample.payload.type = AKU_PAYLOAD_FLOAT;
		sample.payload.float64 = 42;
		
		aku_Status status = aku_write(db, &sample);
		if (status != AKU_SUCCESS && status != AKU_EBUSY) {
			exit(1);
		}
		if (status == AKU_EBUSY) {
			status = aku_write(db, &sample);
			if (status != AKU_SUCCESS) {
				exit(1);
			}
		}
	}
	
	std::cout << "Done writing dummy data\n";
}

int main(int argc, char *argv[]) {
	
	aku_initialize(&panic_handler, &static_logger);
	
	apr_status_t result = aku_create_database("test", "/tmp/db", "/tmp/db", 2, &static_logger);
	if (result != APR_SUCCESS) {
		exit(1);
	}
	
	aku_FineTuneParams params = {
			0, nullptr, 0, 4, 0, 0x1000, 0
	};
	aku_Database* db = aku_open_database("/tmp/db/test.akumuli", params);
	aku_Status status = aku_open_status(db);
	if (status != AKU_SUCCESS) {
		aku_close_database(db);
		exit(1);
	}
	
	writeDummy(db);
	
	aku_close_database(db);
	aku_remove_database("/tmp/test.akumuli", nullptr);
}
