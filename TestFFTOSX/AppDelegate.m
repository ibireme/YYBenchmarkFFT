//
//  AppDelegate.m
//  TestFFTOSX
//
//  Created by ibireme on 14/11/24.
//  Copyright (c) 2014 ibireme. All rights reserved.
//

#import "AppDelegate.h"
#import "FFTRun.h"

@interface AppDelegate ()

@property (weak) IBOutlet NSWindow *window;
@end

@implementation AppDelegate

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification {
    
    dispatch_async(dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0), ^{
        [[FFTRun new] run];
    });
}

- (void)applicationWillTerminate:(NSNotification *)aNotification {
    // Insert code here to tear down your application
}

@end
