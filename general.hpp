//******************** FILE: GENERAL.HPP ********************
//
// Description: This header file contains the definition of functions of general purpose
//
// Author: Vincent Marceau (vincent.marceau.2@ulaval.ca)
//
// Since: October 2012
// Last update: October 2012
//
//*********************************************************

#ifndef GENERAL_HPP_INCLUDED
#define GENERAL_HPP_INCLUDED

// Standard header files
#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>

using namespace std;

typedef unsigned int index_v;

// Function GetDate() definition
// Returns a string with the current date in YYYY-MM-DD format
//
// Input
//     none
// Output
//   - strTime: Current date in YYYY-MM-DD format
//
string GetDate() {
  char cptime[50];
  time_t now = time(NULL);
  strftime(cptime, 50, "%F", localtime(&now)); //uses short month name
  string strTime = cptime;
  return strTime;
} // End function GetDate() definition


// Function GetTime() definition
// Returns a string with the current time in HH:MM:SS format
//
// Input
//     none
// Output
//   - strTime: Current time in HH:MM:SS format
//
string GetTime() {
  char cptime[50];
  time_t now = time(NULL);
  strftime(cptime, 50, "%H:%M:%S", localtime(&now));
  string strTime = cptime;
  return strTime;
} // End function GetTime() definition


// Template function vector_max definition
// Returns the item of maximum value in a container
//
// Input
//     - sequence: container
// Output
//     - maxval: maximum value in sequence
//
template< typename ContentType >
ContentType vector_max(vector<ContentType>& sequence) {
    // Get the maximum element of a given vector

    if (sequence.size() == 0) {
        ContentType maxval = 0;
        return maxval;
    } // end if
    else {
        ContentType maxval = sequence[0];
        for (index_v i=1; i<sequence.size(); i++) {
            if (sequence[i] > maxval)
                maxval = sequence[i];
        } // end for
        return maxval;
    } // end else

} // End template function vector_max definition


// Template function vector_sum definition
// Returns the sum of the items in a container
//
// Input
//     - sequence: container
// Output
//     - sum: sum of sequence's elements
//
template< typename ContentType >
ContentType vector_sum(vector<ContentType>& sequence) {
    // Computes the sum of the elements in a vector

    ContentType sum = 0;
    size_t NumEl = sequence.size();
    for (index_v i=0; i<NumEl; i++)
        sum += sequence[i];
    return sum;
} // End template function vector_sum definition


// Template function vector_mean definition
// Returns the mean value of the items in a container
//
// Input
//     - sequence: container
// Output
//     - sum/NumEl: mean value of the elements in sequence
//
template< typename ContentType >
ContentType vector_mean(vector<ContentType>& sequence) {
    // Computes the mean value of the elements in a given vector

    size_t NumEl = sequence.size();
    ContentType sum = vector_sum(sequence);
    if (NumEl == 0)
        return sum;
    else {
        return sum/NumEl;
    } // end else
} // End template function vector_mean definition


// Template function vector_stddev definition
// Returns the standard deviation of the items in a container
//
// Input
//     - sequence: container
// Output
//     - standard deviation of the elements in sequence
//
template< typename ContentType >
ContentType vector_stddev(vector<ContentType>& sequence) {
    // Computes the standard deviation of the elements in a given vector

    size_t NumEl = sequence.size();
    ContentType sum = 0;
    if (NumEl == 0 || NumEl == 1)
        return sum;
    else {
        ContentType mean = vector_mean(sequence);
        for (index_v i=0; i<NumEl; i++)
            sum = sum + pow((sequence[i]-mean),2.0);
        return sqrt(sum/NumEl);
    } // end else
} // End template function vector_stddev definition


#endif // GENERAL_HPP_INCLUDED

