#include "fesa/model/material.h"

namespace fesa {

Material::Material(int id, char* name) : id_(id), name_(name) {}
Material::~Material() { delete[] name_; }
const ElasticMaterial* Material::AsElastic() const { return nullptr; }
const ThermalMaterial* Material::AsThermal() const { return nullptr; }
int Material::Id() const { return id_; }
const char* Material::Name() const { return name_; }

}  // namespace fesa
