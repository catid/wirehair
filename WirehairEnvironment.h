#pragma once

#include <cstddef>
#include <cstdlib>
#include <cstring>

namespace wirehair {

// Owns an immutable, non-throwing snapshot of an environment value.
class EnvironmentValue
{
public:
    explicit EnvironmentValue(const char* name)
        : Owned(nullptr)
        , Value(nullptr)
    {
#if defined(_MSC_VER)
        std::size_t length = 0;
        if (::_dupenv_s(&Owned, &length, name) == 0) {
            Value = Owned;
        }
#else
        const char* value = std::getenv(name);
        if (value) {
            const std::size_t length = std::strlen(value) + 1;
            Owned = static_cast<char*>(std::malloc(length));
            if (Owned) {
                std::memcpy(Owned, value, length);
                Value = Owned;
            }
        }
#endif
    }

    ~EnvironmentValue()
    {
        std::free(Owned);
    }

    EnvironmentValue(const EnvironmentValue&) = delete;
    EnvironmentValue& operator=(const EnvironmentValue&) = delete;

    bool IsSet() const { return Value != nullptr; }
    const char* Get() const { return Value; }

private:
    char* Owned;
    const char* Value;
};

} // namespace wirehair
