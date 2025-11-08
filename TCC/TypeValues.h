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
inline typename std::enable_if<std::is_floating_point<T>::value, T>::type
falsyValue() {
  return 0.0;
}

// Overload for boolean type
template <>
inline bool falsyValue<bool>() {
  return false;
}

// ------- TRUTHY VALUES ------- 
template <typename T>
inline typename std::enable_if<std::is_integral<T>::value, T>::type
truthyValue() {
  return 1;
}

// Overload for floating-point types 
template <typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, T>::type
truthyValue() {
  return 1.0;
}

// Overload for boolean type
template <>
inline bool truthyValue<bool>() {
  return true;
}

template<typename Test, template<typename...> class Ref>
struct is_specialization : std::false_type {};

template<template<typename...> class Ref, typename... Args>
struct is_specialization<Ref<Args...>, Ref> : std::true_type {};