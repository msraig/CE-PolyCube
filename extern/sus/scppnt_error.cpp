/*! \file scppnt_error.cpp
 \brief Member function definitions for classes used to signal errors.
 
 Definitions of member functions for classes used in throwing
 exceptions indicating various types of run-time errors.
 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

//#ifdef SCPPNT_NO_DIR_PREFIX
#include "scppnt_error.h"
//#else
//#include "scppnt/scppnt_error.h"
//#endif

#include <cstring> // for strlen and strncat
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std
{ using ::strlen; using ::strncat;}
#endif

namespace SCPPNT
{
  Exception::Exception()
  {
    message[0] = 0;
  }

  Exception::~Exception()
  {
  }

  /*
   Concatenate s to the end of the message string
   if there is sufficient space left.
   */
  void Exception::AppendToMessage(const char *s)
  {
    if (s == 0)
      return; // return if pointer is null

    int len = (int)std::strlen(s);
    int currlen = (int)std::strlen(message);
    int freechar = MaxChar() - currlen;

    if (freechar <= 0 || len == 0)
      return;

    if (len > freechar)
      len = freechar;

    std::strncat(message, s, len);

    message[len+currlen] = 0; // add terminating null character
  }

  LogicError::LogicError(const char *type, const char *mess, const char *func)
  {

    /* Create message explaining error */

    AppendToMessage(type);

    if (mess)
    {
      AppendToMessage(mess);
      AppendToMessage("\n");
    }

    if (func)
    {
      AppendToMessage("in function ");
      AppendToMessage(func);
      AppendToMessage("\n");
    }
  }

  LogicError::~LogicError()
  {
  }

  InvalidArgument::InvalidArgument(const char *mess, const char *func) :
    LogicError("Invalid argument: ", mess, func)
  {
  }

  RuntimeError::RuntimeError(const char *mess, const char *func)
  {

    /* Create message explaining error */
    if (mess)
    {
      AppendToMessage(mess);
      AppendToMessage("\n");
    }

    if (func)
    {
      AppendToMessage("in function ");
      AppendToMessage(func);
      AppendToMessage("\n");
    }
  }

  RuntimeError::~RuntimeError()
  {
  }

  BoundsError::BoundsError(const char *func) :
    RuntimeError("Attempted to access element outside valid bounds", func)
  {
  }

  BadDimension::BadDimension(const char *func) :
    RuntimeError("Argument has invalid dimension", func)
  {
  }

} // namespace SCPPNT
