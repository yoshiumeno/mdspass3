#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Input.H>
#include <string>

/**
 * 入力ウィジェットです。
 * ウィジェットの値を登録した外部変数と紐づけ、
 * 外部変数の値の変化がウィジェットに反映されます。
 **/
class CustomInput : public Fl_Input {
public:
    /**
     * 値の種類です。
     * INT_TYPE: 整数
     * FLOAT_TYPE: 浮動小数点数
     * CHAR_TYPE: 文字列
     **/
    enum ValueType {
        INT_TYPE,
        FLOAT_TYPE,
        CHAR_TYPE
    };

    /**
     * コンストラクターです。
     * @param X ウィジェットの x 座標
     * @param Y ウィジェットの y 座標
     * @param W ウィジェットの幅
     * @param H ウィジェットの高さ
     * @param L ウィジェットのラベル
     **/
    CustomInput(int X, int Y, int W, int H, const char* L = 0);
    /* デストラクターです。*/
    virtual ~CustomInput();

    /* 入力と変数を紐づけます。*/
    void attach(int* value);
    void attach(float* value);
    void attach(char* value, size_t size);

    /* 変数の更新を開始します。*/
    void start_polling();
    /* 変数の更新を終了します。*/
    void stop_polling();

private:
    /**
    * 値のポインターです。
    * 整数、浮動小数点数、文字列のいずれかの値を指します。
    **/
    union ValuePtr {
        int* i;
        float* f;
        char* c;
        ValuePtr() : c(nullptr) {}
    };

    /* 値の種類です。*/
    ValueType _value_type;
    /* 値のポインターです。*/
    ValuePtr _value_ptr;
    /* char 配列のサイズです。*/
    size_t _char_array_size;

    /* 変数の更新を行うコールバック関数です。*/
    static void polling_callback(void* v);

    /**
     * イベント ハンドラーです。
     * @param event イベント識別子
     **/
    int handle(int event);
    /* 入力値から変数を更新します。*/
    void update_variable_from_input();
};
