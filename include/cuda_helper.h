#include "stdio.h"
#ifndef _CUDA_HELPER_H_
#define _CUDA_HELPER_H_

#define checkCudaErrors(val)           check ( (val), #val, __FILE__, __LINE__ )

#ifndef MIN
#define MIN(a,b) ((a < b) ? a : b)
#endif
#ifndef MAX
#define MAX(a,b) (a > b ? a : b)
#endif

#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#ifndef STRCASECMP
#define STRCASECMP  _stricmp
#endif
#ifndef STRNCASECMP
#define STRNCASECMP _strnicmp
#endif
#ifndef STRCPY
#define STRCPY(sFilePath, nLength, sPath) strcpy_s(sFilePath, nLength, sPath)
#endif

#ifndef FOPEN
#define FOPEN(fHandle,filename,mode) fopen_s(&fHandle, filename, mode)
#endif
#ifndef FOPEN_FAIL
#define FOPEN_FAIL(result) (result != 0)
#endif
#ifndef SSCANF
#define SSCANF sscanf_s
#endif
#ifndef SPRINTF
#define SPRINTF sprintf_s
#endif
#else // Linux Includes
#include <string.h>
#include <strings.h>

#ifndef STRCASECMP
#define STRCASECMP  strcasecmp
#endif
#ifndef STRNCASECMP
#define STRNCASECMP strncasecmp
#endif
#ifndef STRCPY
#define STRCPY(sFilePath, nLength, sPath) strcpy(sFilePath, sPath)
#endif

#ifndef FOPEN
#define FOPEN(fHandle,filename,mode) (fHandle = fopen(filename, mode))
#endif
#ifndef FOPEN_FAIL
#define FOPEN_FAIL(result) (result == NULL)
#endif
#ifndef SSCANF
#define SSCANF sscanf
#endif
#ifndef SPRINTF
#define SPRINTF sprintf
#endif
#endif

#ifndef EXIT_WAIVED
#define EXIT_WAIVED 2
#endif

#include <cuda_runtime_api.h>
#include <cuda.h>


template< typename T >
void check(T result, char const *const func, const char *const file, int const line)
{
	if (result)
	{
		printf( "CUDA error at %s:%d code=%d(%s) \"%s\" \n",
			file, line, static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);
		cudaDeviceReset();
		// Make sure we call CUDA Device Reset before exiting
		exit(EXIT_FAILURE);
	}
}
#ifdef __DRIVER_TYPES_H__
static const char *_cudaGetErrorEnum(cudaError_t error)
{
	switch (error)
	{
	case cudaSuccess:
		return "cudaSuccess";

	case cudaErrorMissingConfiguration:
		return "cudaErrorMissingConfiguration";

	case cudaErrorMemoryAllocation:
		return "cudaErrorMemoryAllocation";

	case cudaErrorInitializationError:
		return "cudaErrorInitializationError";

	case cudaErrorLaunchFailure:
		return "cudaErrorLaunchFailure";

	case cudaErrorPriorLaunchFailure:
		return "cudaErrorPriorLaunchFailure";

	case cudaErrorLaunchTimeout:
		return "cudaErrorLaunchTimeout";

	case cudaErrorLaunchOutOfResources:
		return "cudaErrorLaunchOutOfResources";

	case cudaErrorInvalidDeviceFunction:
		return "cudaErrorInvalidDeviceFunction";

	case cudaErrorInvalidConfiguration:
		return "cudaErrorInvalidConfiguration";

	case cudaErrorInvalidDevice:
		return "cudaErrorInvalidDevice";

	case cudaErrorInvalidValue:
		return "cudaErrorInvalidValue";

	case cudaErrorInvalidPitchValue:
		return "cudaErrorInvalidPitchValue";

	case cudaErrorInvalidSymbol:
		return "cudaErrorInvalidSymbol";

	case cudaErrorMapBufferObjectFailed:
		return "cudaErrorMapBufferObjectFailed";

	case cudaErrorUnmapBufferObjectFailed:
		return "cudaErrorUnmapBufferObjectFailed";

	case cudaErrorInvalidHostPointer:
		return "cudaErrorInvalidHostPointer";

	case cudaErrorInvalidDevicePointer:
		return "cudaErrorInvalidDevicePointer";

	case cudaErrorInvalidTexture:
		return "cudaErrorInvalidTexture";

	case cudaErrorInvalidTextureBinding:
		return "cudaErrorInvalidTextureBinding";

	case cudaErrorInvalidChannelDescriptor:
		return "cudaErrorInvalidChannelDescriptor";

	case cudaErrorInvalidMemcpyDirection:
		return "cudaErrorInvalidMemcpyDirection";

	case cudaErrorAddressOfConstant:
		return "cudaErrorAddressOfConstant";

	case cudaErrorTextureFetchFailed:
		return "cudaErrorTextureFetchFailed";

	case cudaErrorTextureNotBound:
		return "cudaErrorTextureNotBound";

	case cudaErrorSynchronizationError:
		return "cudaErrorSynchronizationError";

	case cudaErrorInvalidFilterSetting:
		return "cudaErrorInvalidFilterSetting";

	case cudaErrorInvalidNormSetting:
		return "cudaErrorInvalidNormSetting";

	case cudaErrorMixedDeviceExecution:
		return "cudaErrorMixedDeviceExecution";

	case cudaErrorCudartUnloading:
		return "cudaErrorCudartUnloading";

	case cudaErrorUnknown:
		return "cudaErrorUnknown";

	case cudaErrorNotYetImplemented:
		return "cudaErrorNotYetImplemented";

	case cudaErrorMemoryValueTooLarge:
		return "cudaErrorMemoryValueTooLarge";

	case cudaErrorInvalidResourceHandle:
		return "cudaErrorInvalidResourceHandle";

	case cudaErrorNotReady:
		return "cudaErrorNotReady";

	case cudaErrorInsufficientDriver:
		return "cudaErrorInsufficientDriver";

	case cudaErrorSetOnActiveProcess:
		return "cudaErrorSetOnActiveProcess";

	case cudaErrorInvalidSurface:
		return "cudaErrorInvalidSurface";

	case cudaErrorNoDevice:
		return "cudaErrorNoDevice";

	case cudaErrorECCUncorrectable:
		return "cudaErrorECCUncorrectable";

	case cudaErrorSharedObjectSymbolNotFound:
		return "cudaErrorSharedObjectSymbolNotFound";

	case cudaErrorSharedObjectInitFailed:
		return "cudaErrorSharedObjectInitFailed";

	case cudaErrorUnsupportedLimit:
		return "cudaErrorUnsupportedLimit";

	case cudaErrorDuplicateVariableName:
		return "cudaErrorDuplicateVariableName";

	case cudaErrorDuplicateTextureName:
		return "cudaErrorDuplicateTextureName";

	case cudaErrorDuplicateSurfaceName:
		return "cudaErrorDuplicateSurfaceName";

	case cudaErrorDevicesUnavailable:
		return "cudaErrorDevicesUnavailable";

	case cudaErrorInvalidKernelImage:
		return "cudaErrorInvalidKernelImage";

	case cudaErrorNoKernelImageForDevice:
		return "cudaErrorNoKernelImageForDevice";

	case cudaErrorIncompatibleDriverContext:
		return "cudaErrorIncompatibleDriverContext";

	case cudaErrorPeerAccessAlreadyEnabled:
		return "cudaErrorPeerAccessAlreadyEnabled";

	case cudaErrorPeerAccessNotEnabled:
		return "cudaErrorPeerAccessNotEnabled";

	case cudaErrorDeviceAlreadyInUse:
		return "cudaErrorDeviceAlreadyInUse";

	case cudaErrorProfilerDisabled:
		return "cudaErrorProfilerDisabled";

	case cudaErrorProfilerNotInitialized:
		return "cudaErrorProfilerNotInitialized";

	case cudaErrorProfilerAlreadyStarted:
		return "cudaErrorProfilerAlreadyStarted";

	case cudaErrorProfilerAlreadyStopped:
		return "cudaErrorProfilerAlreadyStopped";

		/* Since CUDA 4.0*/
	case cudaErrorAssert:
		return "cudaErrorAssert";

	case cudaErrorTooManyPeers:
		return "cudaErrorTooManyPeers";

	case cudaErrorHostMemoryAlreadyRegistered:
		return "cudaErrorHostMemoryAlreadyRegistered";

	case cudaErrorHostMemoryNotRegistered:
		return "cudaErrorHostMemoryNotRegistered";

		/* Since CUDA 5.0 */
	case cudaErrorOperatingSystem:
		return "cudaErrorOperatingSystem";

		//        case cudaErrorPeerAccssUnsupported:
		//            return "cudaErrorPeerAccessUnsupported";

	case cudaErrorLaunchMaxDepthExceeded:
		return "cudaErrorLaunchMaxDepthExceeded";

	case cudaErrorLaunchFileScopedTex:
		return "cudaErrorLaunchFileScopedTex";

	case cudaErrorLaunchFileScopedSurf:
		return "cudaErrorLaunchFileScopedSurf";

	case cudaErrorSyncDepthExceeded:
		return "cudaErrorSyncDepthExceeded";

	case cudaErrorLaunchPendingCountExceeded:
		return "cudaErrorLaunchPendingCountExceeded";

	case cudaErrorNotPermitted:
		return "cudaErrorNotPermitted";

	case cudaErrorNotSupported:
		return "cudaErrorNotSupported";

		/* Since CUDA 6.0 */
	case cudaErrorHardwareStackError:
		return "cudaErrorHardwareStackError";

	case cudaErrorIllegalInstruction:
		return "cudaErrorIllegalInstruction";

	case cudaErrorMisalignedAddress:
		return "cudaErrorMisalignedAddress";

	case cudaErrorInvalidAddressSpace:
		return "cudaErrorInvalidAddressSpace";

	case cudaErrorInvalidPc:
		return "cudaErrorInvalidPc";

	case cudaErrorIllegalAddress:
		return "cudaErrorIllegalAddress";

		/* Since CUDA 6.5*/
	case cudaErrorInvalidPtx:
		return "cudaErrorInvalidPtx";

	case cudaErrorInvalidGraphicsContext:
		return "cudaErrorInvalidGraphicsContext";

	case cudaErrorStartupFailure:
		return "cudaErrorStartupFailure";

	case cudaErrorApiFailureBase:
		return "cudaErrorApiFailureBase";

		/* Since CUDA 8.0*/
	case cudaErrorNvlinkUncorrectable:
		return "cudaErrorNvlinkUncorrectable";
	}

	return "<unknown>";
}

inline int stringRemoveDelimiter(char delimiter, const char *string)
{
	int string_start = 0;

	while (string[string_start] == delimiter)
	{
		string_start++;
	}

	if (string_start >= (int)strlen(string) - 1)
	{
		return 0;
	}

	return string_start;
}

inline int getFileExtension(char *filename, char **extension)
{
	int string_length = (int)strlen(filename);

	while (filename[string_length--] != '.')
	{
		if (string_length == 0)
			break;
	}

	if (string_length > 0) string_length += 2;

	if (string_length == 0)
		*extension = NULL;
	else
		*extension = &filename[string_length];

	return string_length;
}
inline bool checkCmdLineFlag(const int argc, const char **argv, const char *string_ref)
{
	bool bFound = false;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			const char *string_argv = &argv[i][string_start];

			const char *equal_pos = strchr(string_argv, '=');
			int argv_length = (int)(equal_pos == 0 ? strlen(string_argv) : equal_pos - string_argv);

			int length = (int)strlen(string_ref);

			if (length == argv_length && !STRNCASECMP(string_argv, string_ref, length))
			{
				bFound = true;
				continue;
			}
		}
	}

	return bFound;
}

// This function wraps the CUDA Driver API into a template function
template <class T>
inline bool getCmdLineArgumentValue(const int argc, const char **argv, const char *string_ref, T *value)
{
	bool bFound = false;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			const char *string_argv = &argv[i][string_start];
			int length = (int)strlen(string_ref);

			if (!STRNCASECMP(string_argv, string_ref, length))
			{
				if (length + 1 <= (int)strlen(string_argv))
				{
					int auto_inc = (string_argv[length] == '=') ? 1 : 0;
					*value = (T)atoi(&string_argv[length + auto_inc]);
				}

				bFound = true;
				i = argc;
			}
		}
	}

	return bFound;
}

inline int getCmdLineArgumentInt(const int argc, const char **argv, const char *string_ref)
{
	bool bFound = false;
	int value = -1;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			const char *string_argv = &argv[i][string_start];
			int length = (int)strlen(string_ref);

			if (!STRNCASECMP(string_argv, string_ref, length))
			{
				if (length + 1 <= (int)strlen(string_argv))
				{
					int auto_inc = (string_argv[length] == '=') ? 1 : 0;
					value = atoi(&string_argv[length + auto_inc]);
				}
				else
				{
					value = 0;
				}

				bFound = true;
				continue;
			}
		}
	}

	if (bFound)
	{
		return value;
	}
	else
	{
		return 0;
	}
}

inline float getCmdLineArgumentFloat(const int argc, const char **argv, const char *string_ref)
{
	bool bFound = false;
	float value = -1;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			const char *string_argv = &argv[i][string_start];
			int length = (int)strlen(string_ref);

			if (!STRNCASECMP(string_argv, string_ref, length))
			{
				if (length + 1 <= (int)strlen(string_argv))
				{
					int auto_inc = (string_argv[length] == '=') ? 1 : 0;
					value = (float)atof(&string_argv[length + auto_inc]);
				}
				else
				{
					value = 0.f;
				}

				bFound = true;
				continue;
			}
		}
	}

	if (bFound)
	{
		return value;
	}
	else
	{
		return 0;
	}
}

inline bool getCmdLineArgumentString(const int argc, const char **argv,
	const char *string_ref, char **string_retval)
{
	bool bFound = false;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			char *string_argv = (char *)&argv[i][string_start];
			int length = (int)strlen(string_ref);

			if (!STRNCASECMP(string_argv, string_ref, length))
			{
				*string_retval = &string_argv[length + 1];
				bFound = true;
				continue;
			}
		}
	}

	if (!bFound)
	{
		*string_retval = NULL;
	}

	return bFound;
}
#endif
#endif
