#include <Arduino.h>
#include <arduinoFFT.h>  
#include <math.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include <algorithm>

const int ADC_GPIO = 5;   // GPIO5 для ADC1_CH4
const int SDA_PIN = 17;
const int SCL_PIN = 18;

const uint32_t SAMPLE_RATE = 10000; //для интрепр 5000 но тогда мы знаем что поедет
const uint16_t FFT_SIZE = 2048; //для интерполяции может надо 4096 но поедут другие штуки
const uint8_t MAX_HARMONIC = 3; //поизменять глянуть
const uint16_t RMS_BLOCK = 512;
const double RMS_DB_THRESHOLD = -25.0;
const int CALIBRATE_SAMPLES = 512;

#define SCREEN_WIDTH 128
#define SCREEN_HEIGHT 64
Adafruit_SSD1306 display(SCREEN_WIDTH, SCREEN_HEIGHT, &Wire, -1);

ArduinoFFT<double>* FFT = nullptr;
double *vReal = nullptr;
double *vImag = nullptr;
double *spec = nullptr;
double *hps = nullptr;
double adcCenter = 2048.0;

const uint8_t MEDIAN_SIZE = 7;
double freqHistory[MEDIAN_SIZE] = {0};
uint8_t histIndex = 0;



struct GuitarString {
  const char* name;     
  double freq;         
  double lowBound;    
  double highBound;   
};

const GuitarString CLASSIC_TUNING[] = {
  {"1 (E)", 329.63, 320.0, 340.0}, 
  {"2 (B)", 246.94, 240.0, 254.0},
  {"3 (G)", 196.00, 190.0, 202.0},
  {"4 (D)", 146.83, 142.0, 152.0},
  {"5 (A)", 110.00, 106.0, 114.0},
  {"6 (E)", 82.41,  80.0,  85.0}    
};

const int STRING_COUNT = 6;
const double TUNING_TOLERANCE = 0.5; 

double parabolic_interpolation(const double *mag, int k, int Nhalf);
double to_dbfs(double rms, double fullscale);
double medianFilter(double newFreq);
int detectString(double frequency);
String getTuningAdvice(double currentFreq, double targetFreq);

void setup() {
  Serial.begin(115200);
  delay(500);
  
  Wire.begin(SDA_PIN, SCL_PIN);
  if (!display.begin(SSD1306_SWITCHCAPVCC, 0x3C)) {
    Serial.println("SSD1306 allocation failed");
    while (1) delay(10);
  }
  
  display.clearDisplay();
  display.setTextColor(WHITE);
  display.setTextSize(1);
  display.setCursor(0, 0);
  display.println("Starting...");
  display.display();
  

  vReal = (double*)malloc(sizeof(double) * FFT_SIZE);
  vImag = (double*)malloc(sizeof(double) * FFT_SIZE);
  spec = (double*)malloc(sizeof(double) * (FFT_SIZE/2));
  hps = (double*)malloc(sizeof(double) * (FFT_SIZE/2));
  
  if (!vReal || !vImag || !spec || !hps) {
    Serial.println("Memory allocation failed!");
    while(1);
  }
  
  Serial.print("Free heap: ");
  Serial.println(ESP.getFreeHeap());
  

  analogReadResolution(12);         
  analogSetAttenuation(ADC_11db);  
  

  Serial.println("Calibrating ADC center...");
  display.clearDisplay();
  display.setCursor(0, 0);
  display.println("Calibrating...");
  display.display();
  
  double sum = 0;
  uint32_t calStartTime = micros();
  for (int i = 0; i < CALIBRATE_SAMPLES; i++) {
    sum += analogRead(ADC_GPIO);

    while (micros() - calStartTime < i * (1000000.0 / SAMPLE_RATE)) {}
  }
  adcCenter = sum / CALIBRATE_SAMPLES;
  
  Serial.print("ADC center: ");
  Serial.println(adcCenter, 3);
  

  Serial.println("Testing ADC...");
  for (int i = 0; i < 10; i++) {
    Serial.print("ADC reading ");
    Serial.print(i);
    Serial.print(": ");
    Serial.println(analogRead(ADC_GPIO));
    delay(100);
  }
  
  FFT = new ArduinoFFT<double>(vReal, vImag, FFT_SIZE, SAMPLE_RATE);
  if (FFT == nullptr) {
    Serial.println("Failed to create FFT object");
    while(1);
  }
  
  Serial.println("Tuner ready");
  display.clearDisplay();
  display.setCursor(0, 0);
  display.println("Classic Tuning");
  display.println("Ready...");
  display.display();
  delay(2000);
}

void loop() {
  double sum = 0;
  uint32_t rmsStartTime = micros();
  
  for (uint16_t i = 0; i < RMS_BLOCK; i++) {
    double s = analogRead(ADC_GPIO) - adcCenter;
    sum += s * s;
    while (micros() - rmsStartTime < i * (1000000.0 / SAMPLE_RATE)) {}
  }
  
  double rms = sqrt(sum / RMS_BLOCK);
  double dbfs = to_dbfs(rms, adcCenter);
  
  if (dbfs < RMS_DB_THRESHOLD) {
    display.clearDisplay();
    display.setCursor(0, 0);
    display.println("No signal");
    display.display();
    delay(100);
    return;
  }
  

  uint32_t sampleStartTime = micros();
  
  for (uint16_t i = 0; i < FFT_SIZE; i++) {
    vReal[i] = analogRead(ADC_GPIO) - adcCenter;
    vImag[i] = 0.0;
    while (micros() - sampleStartTime < i * (1000000.0 / SAMPLE_RATE)) {}
  }
  

  uint32_t actualDuration = micros() - sampleStartTime;
  float actualFs = (FFT_SIZE * 1000000.0) / actualDuration;
  


  FFT->windowing(FFTWindow::Hamming, FFTDirection::Forward);
  FFT->compute(FFTDirection::Forward);
  FFT->complexToMagnitude();
  vReal[0] = 0.0;  
  

  uint32_t half = FFT_SIZE / 2;
  double deltaF = (double)SAMPLE_RATE / FFT_SIZE;
  

  double maxSpec = 1e-12;
  for (uint32_t i = 0; i < half; i++) {
    spec[i] = vReal[i] + 1e-12;
    if (spec[i] > maxSpec) maxSpec = spec[i];
  }
  
  for (uint32_t i = 0; i < half; i++) {
    spec[i] /= maxSpec;
  }
  

  for (uint32_t i = 0; i < half; i++) {
    hps[i] = log(spec[i] + 1e-12);
  }
  
  for (uint8_t k = 2; k <= MAX_HARMONIC; k++) {
    for (uint32_t i = 0; i < half / k; i++) {
      hps[i] += log(spec[i * k] + 1e-12);
    }
  }
  

  uint32_t minIdx = max(1, (int)round(70.0 / deltaF));   
  uint32_t maxIdx = min((int)half - 1, (int)round(400.0 / deltaF));
  
  double maxHps = -1e300;
  int maxHpsIdx = minIdx;
  
  for (uint32_t i = minIdx; i <= maxIdx; i++) {
    if (hps[i] > maxHps) {
      maxHps = hps[i];
      maxHpsIdx = i;
    }
  }

  double rawFreq = maxHpsIdx * deltaF;
  

  double fineFreq = parabolic_interpolation(hps, maxHpsIdx, half) * deltaF;
  double stableFreq = medianFilter(fineFreq);
  
  int detectedString = detectString(stableFreq);
  String tuningAdvice = "";
  double freqDiff = 0;
  
  if (detectedString >= 0) {
    freqDiff = stableFreq - CLASSIC_TUNING[detectedString].freq;
    tuningAdvice = getTuningAdvice(stableFreq, CLASSIC_TUNING[detectedString].freq);
  }
  
  Serial.print("Freq: ");
  Serial.print(stableFreq, 1);
  Serial.print(" Hz | String: ");
  if (detectedString >= 0) {
    Serial.print(CLASSIC_TUNING[detectedString].name);
  } else {
    Serial.print("Unknown");
  }
  Serial.print(" | Diff: ");
  Serial.print(freqDiff, 2);
  Serial.print(" Hz | Advice: ");
  Serial.println(tuningAdvice);
  
  display.clearDisplay();
  
  display.setTextSize(2);
  display.setCursor(0, 0);
  display.print(stableFreq, 1);
  display.print("Hz");
  
  if (detectedString >= 0) {
    String stringNum = String(CLASSIC_TUNING[detectedString].name);
    int spacePos = stringNum.indexOf(' ');
    if (spacePos > 0) {
      stringNum = stringNum.substring(0, spacePos);
    }
    
    display.setTextSize(2); 
    display.setCursor(90, 0); 
    display.print(stringNum);
  }
  // display.setTextSize(1);
  // display.setCursor(0, 25);
  // if (detectedString >= 0) {
  //   display.print("String ");
  //   display.print(CLASSIC_TUNING[detectedString].name);
  // } else {
  //   display.print("String: Unknown");
  // }
  
  display.setTextSize(2);
  display.setCursor(0, 30);
  display.print("Diff:");
  if (freqDiff > 0) {
    display.print("+");
  }
  display.print(freqDiff, 1);
  
  display.setTextSize(1);
  display.setCursor(0, 50);
  if (detectedString >= 0) {
    if (fabs(freqDiff) < TUNING_TOLERANCE) {
      display.print("PERFECT! ✓");
    } else {
      display.print(tuningAdvice);
    }
  }
  
  // display.setCursor(0, 55);
  // display.print("dBFS: ");
  // display.print(dbfs, 1);
  
  display.display();
  
  delay(2000);  
}

int detectString(double frequency) {
  for (int i = 0; i < STRING_COUNT; i++) {
    if (frequency >= CLASSIC_TUNING[i].lowBound && 
        frequency <= CLASSIC_TUNING[i].highBound) {
      return i;
    }
  }
  
  int closestString = 0;
  double minDiff = 1000.0;
  
  for (int i = 0; i < STRING_COUNT; i++) {
    double diff = fabs(frequency - CLASSIC_TUNING[i].freq);
    if (diff < minDiff) {
      minDiff = diff;
      closestString = i;
    }
  }
  
  return closestString;
}

// Функция для получения рекомендации по настройке
String getTuningAdvice(double currentFreq, double targetFreq) {
  double diff = currentFreq - targetFreq;
  double absDiff = fabs(diff);
  
  if (absDiff < 0.2) {
    return "Perfect!";
  } else if (diff > 0) {
    if (absDiff > 5.0) return "Loosen a lot";
    else if (absDiff > 2.0) return "Loosen more";
    else return "Loosen slightly";
  } else {
    if (absDiff > 5.0) return "Tighten a lot";
    else if (absDiff > 2.0) return "Tighten more";
    else return "Tighten slightly";
  }
}

double medianFilter(double newFreq) {
  static int initialized = 0;
  static double lastValidFreq = 0;
  
  // скользящее среднее
  if (initialized < MEDIAN_SIZE) {
    freqHistory[initialized] = newFreq;
    initialized++;
    
    // Среднее арифметическое для начала записи
    double sum = 0;
    for(int i = 0; i < initialized; i++) {
      sum += freqHistory[i];
    }
    lastValidFreq = sum / initialized;
    
    if (initialized == MEDIAN_SIZE) {
      Serial.print("Filter initialized. Initial median: ");
      Serial.println(lastValidFreq, 1);
    }
    
    return lastValidFreq;
  }
  
  // Проверка на смену струны
  if (lastValidFreq > 0 && fabs(newFreq - lastValidFreq) > 20.0) {
    for(int i = 0; i < MEDIAN_SIZE; i++) {
      freqHistory[i] = newFreq;
    }
    histIndex = 0;
    
    double medianAfterReset = newFreq;  // сначала все значения будут одинаковые
    
    Serial.print("STRING CHANGED: ");
    Serial.print(lastValidFreq, 1);
    Serial.print(" Hz -> ");
    Serial.print(newFreq, 1);
    Serial.print(" Hz | Median: ");
    Serial.println(medianAfterReset, 1);
    
    lastValidFreq = medianAfterReset;
    return medianAfterReset;
  }
  
  freqHistory[histIndex] = newFreq;
  histIndex = (histIndex + 1) % MEDIAN_SIZE;
  
  double temp[MEDIAN_SIZE];
  memcpy(temp, freqHistory, sizeof(temp));
  
  for (int i = 0; i < MEDIAN_SIZE - 1; i++) {
    for (int j = i + 1; j < MEDIAN_SIZE; j++) {
      if (temp[i] > temp[j]) {
        double swap = temp[i];
        temp[i] = temp[j];
        temp[j] = swap;
      }
    }
  }
  
  double currentMedian = temp[MEDIAN_SIZE / 2];
  
  if (fabs(currentMedian - newFreq) > 10.0) {
    Serial.print("Filter correction: ");
    Serial.print(newFreq, 1);
    Serial.print(" Hz -> Median: ");
    Serial.println(currentMedian, 1);
  }
  
  lastValidFreq = currentMedian;
  return currentMedian;
}

double parabolic_interpolation(const double *mag, int k, int Nhalf) {
  if (k <= 0 || k >= Nhalf - 1) return (double)k;
  
  double alpha = mag[k - 1];
  double beta = mag[k];
  double gamma = mag[k + 1];
  double denom = alpha - 2.0 * beta + gamma;
  
  if (fabs(denom) < 1e-12) return (double)k;
  return k + 0.5 * (alpha - gamma) / denom;
}

double to_dbfs(double rms, double fullscale) {
  return 20.0 * log10((rms + 1e-12) / fullscale);
}

void cleanup() {
  if (FFT != nullptr) {
    delete FFT;
    FFT = nullptr;
  }
  
  if (vReal != nullptr) free(vReal);
  if (vImag != nullptr) free(vImag);
  if (spec != nullptr) free(spec);
  if (hps != nullptr) free(hps);
  
  Serial.println("Memory cleaned up");
}
