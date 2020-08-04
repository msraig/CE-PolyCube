/*! \file scppnt_error.h
 \brief Exception classes used to signal errors.
 
 Declarations of classes used in throwing
 exceptions indicating various types of run-time errors.
 */

/*

 Simple C++ Numerical Toolkit (SCPPNT)
 http://www.smallwaters.com/software/cpp/scppnt.html
 This release updates original work contributed by 
 Brad Hanson (http://www.b-a-h.com/) in 2001.

 */

#ifndef SCPPNTERROR
#define SCPPNTERROR

namespace SCPPNT
{

  //! General exception class from which classes indicating more specific types are errors are derived.
  class Exception
  {

public:

    //! Constructor
    Exception();

    //! Destructor
    virtual ~Exception();

    //! Returns message describing the error that occurred
    virtual const char *what() const
    {
      return message;
    }

protected:

    //! concatenate string to end of message if there is sufficient space left
    void AppendToMessage(const char *s);

    //! return maximum number of characters in message
    int MaxChar()
    {
      return 200;
    }

    //! space to hold message
    char message[201];

  };

  //! Errors in program logic that should be caught during program debugging.
  class LogicError : public Exception
  {

public:

    //! Constructor from error type, error message, and function name form which exception is thrown
    LogicError(const char *type, const char *message, const char *func);

    //! Destructor
    virtual ~LogicError();

  };

  //! A function argument is invalid
  class InvalidArgument : public LogicError
  {
public:

    //! Constructor from error message and function name form which exception is thrown
    InvalidArgument(const char *message, const char *func);

  };

  //! General run-time errors. Classes for specific types of run-time errors are derived from this class.
  class RuntimeError : public Exception
  {
public:

    //! Constructor from error message and function name form which exception is thrown
    RuntimeError(const char *message, const char *func);

    //! Destructor
    virtual ~RuntimeError();
  };

  //! Attempt to access element outside vector or matrix bounds
  class BoundsError : public RuntimeError
  {
public:

    //! Constructor from function name from which exception is thrown
    BoundsError(const char *func);

  };

  //! Bad dimension of vector or matrix argument of function
  class BadDimension : public RuntimeError
  {
public:

    //! Constructor from function name from which exception is thrown
    BadDimension(const char *func);

  };

} // namespace SCPPNT

#endif // SCPPNTERROR
