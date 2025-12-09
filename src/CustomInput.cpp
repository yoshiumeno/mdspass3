#include "CustomInput.h"
#include <sstream>
#include <iomanip>
#include "myheader.h"

CustomInput::CustomInput(int X, int Y, int W, int H, const char* L)
    : Fl_Input(X, Y, W, H, L), _value_type(CHAR_TYPE), _char_array_size(0)
{
    // 値のポインターを初期化します。
    _value_ptr.c = nullptr;
}

CustomInput::~CustomInput()
{
    // 値の更新を終了します。
    stop_polling();
}

void CustomInput::attach(int* value)
{
    // 整数の値を紐づけます。
    _value_type = INT_TYPE;
    _value_ptr.i = value;
}

void CustomInput::attach(float* value)
{
    // 浮動小数点数の値を紐づけます。
    _value_type = FLOAT_TYPE;
    _value_ptr.f = value;
}

void CustomInput::attach(char* value, size_t size)
{
    // 文字列の値を紐づけます。
    _value_type = CHAR_TYPE;
    _value_ptr.c = value;
    _char_array_size = size;
}

void CustomInput::start_polling()
{
    // 値の更新を開始します。
    Fl::add_timeout(POLLING_INTERVAL, polling_callback, this);
}

void CustomInput::stop_polling()
{
    // 値の更新を終了します。
    Fl::remove_timeout(polling_callback, this);
}

void CustomInput::polling_callback(void* v)
{
    // 引数を元の型にキャストします。
    CustomInput* widget = static_cast<CustomInput*>( v );
    // 文字列ストリームを作成します。
    std::ostringstream string_stream;

    // 外部変数の値を文字列に変換します。
    // 値の種類に応じて、処理を行います。
    switch (widget->_value_type)
    {
        // 整数の場合です。
        case INT_TYPE:
            // 整数値を文字列に変換します。
            if (widget->_value_ptr.i) string_stream << *( widget->_value_ptr.i );
            break;
            // 浮動小数点数の場合です。
        case FLOAT_TYPE:
            // 小数点以下 3 桁までを文字列にして表示します。
            if (widget->_value_ptr.f) string_stream << std::fixed << std::setprecision(3) << *( widget->_value_ptr.f );
            break;
            // 文字列の場合です。
        case CHAR_TYPE:
            if (widget->_value_ptr.c) string_stream << std::string(widget->_value_ptr.c);
            break;
    }
    // 外部変数の値が、入力フィールドと異なる場合です。
    if (string_stream.str() != widget->value())
    {
        // 入力フィールドの値を更新します。
        widget->value(string_stream.str().c_str());
    }
    // 次の更新をスケジュールします。
    Fl::repeat_timeout(POLLING_INTERVAL, polling_callback, widget);
}

int CustomInput::handle(int event)
{
    switch (event)
    {
        case FL_FOCUS:
            // 入力開始時 (フォーカス取得時) に自動更新を停止します。
            stop_polling();
            break;
        case FL_UNFOCUS:
            // 入力終了時 (フォーカス失った時) に変数を更新し、自動更新を再開します。
            update_variable_from_input();
            start_polling();
            break;
        default:
            break;
    }
    // FLTK のデフォルトのイベント ハンドリングを継続します。
    return Fl_Input::handle(event);
}

void CustomInput::update_variable_from_input()
{
    std::string input_value = this->value();
    switch (_value_type)
    {
        case INT_TYPE:
            try
            {
                int new_value = std::stoi(input_value);
                // 有効な入力の場合、値を更新します。
                *_value_ptr.i = new_value;
            }
            catch (const std::invalid_argument&)
            {
                // 無効な入力の場合、元の値を表示します。
                this->value(std::to_string(*_value_ptr.i).c_str());
            }
            catch (const std::out_of_range&)
            {
                // 入力値が範囲外の場合も、元の値を表示します。
                this->value(std::to_string(*_value_ptr.i).c_str());
            }
            break;
        case FLOAT_TYPE:
            try
            {
                //float new_value = std::stoi(input_value);
                float new_value = std::stof(input_value); //YU2024.08
                // 有効な入力の場合、値を更新します。
                *_value_ptr.f = new_value;
            }
            catch (const std::invalid_argument&)
            {
                // 無効な入力の場合、元の値を表示します。
                this->value(std::to_string(*_value_ptr.f).c_str());
            }
            catch (const std::out_of_range&)
            {
                // 入力値が範囲外の場合も、元の値を表示します。
                this->value(std::to_string(*_value_ptr.f).c_str());
            }
            break;
        case CHAR_TYPE:
            if (_value_ptr.c && _char_array_size > 0)
            {
                // input_value の長さが配列のサイズを超えないように制限します。
                size_t copyLength = std::min(input_value.length(), _char_array_size - 1);
                // 必要な分だけコピーし、ヌル終端を保証します。
                strncpy(_value_ptr.c, input_value.c_str(), copyLength);
                _value_ptr.c[copyLength] = '\0';
            }
            break;
    }
}
