#pragma once
#include <type_traits>

// ------- FALSY VALUES -------  
template <typename T>
typename std::enable_if<std::is_integral<T>::value, T>::type 
falsyValue() {
  static_assert(!std::is_same<T, bool>::value, "only numbers!");
  return 0;
}

// Overload for floating-point types 
template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
falsyValue() {
  return 0.0;
}

// Overload for boolean type
template <>
bool falsyValue<bool>() {
  return false;
}

// ------- TRUTHY VALUES ------- 
template <typename T>
typename std::enable_if<std::is_integral<T>::value, T>::type
truthyValue() {
  return 1;
}

// Overload for floating-point types 
template <typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
truthyValue() {
  return 1.0;
}

// Overload for boolean type
template <>
bool truthyValue<bool>() {
  return true;
}