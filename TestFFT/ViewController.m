//
//  ViewController.m
//  TestFFT
//
//  Created by ibireme on 14/11/24.
//  Copyright (c) 2014 ibireme. All rights reserved.
//

#import "ViewController.h"
#import "FFTRun.h"

@interface ViewController ()

@end

@implementation ViewController

- (void)viewDidLoad {
    [super viewDidLoad];
    
    
    dispatch_after(dispatch_time(DISPATCH_TIME_NOW, (int64_t)(1.0 * NSEC_PER_SEC)), dispatch_get_main_queue(), ^{
        dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^{
            
            [[FFTRun new] run];
            
        });
    });
    
}



@end
